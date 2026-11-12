#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "lndcal.h"
#include "keyvalue.h"
#include "const.h"
#include "param.h"
#include "input.h"
#include "lut.h"
#include "output.h"
#include "cal.h"
#include "error.h"

#include <time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

int main (int argc, char *argv[]) {
  Param_t *param = NULL;
  Input_t *input = NULL;
  Lut_t *lut = NULL;
  Output_t *output = NULL;
  Output_t *output_th = NULL;
  int iline, ib;
  int curr_line;        /* current location in the QA band for the start
                           of the current line and for the current pixel */
  unsigned char *line_in = NULL;
  int16 *line_in_sun_zen = NULL;     /* solar zenith representative band */
  uint16_t *line_qa = NULL;
  uint16_t *line_out = NULL;
  uint16_t *line_out_th = NULL;
  int nps,nls, nps6, nls6;
  int odometer_flag=0;
  char msgbuf[1024];
  char envi_file[STR_SIZE]; /* name of the output ENVI header file */
  char *cptr=NULL;          /* pointer to the file extension */
  size_t input_psize;
  int nband_refl = NBAND_REFL_MAX;
  Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
  Envi_header_t envi_hdr;   /* output ENVI header information */

  /* Read the parameters from the input parameter file */
  param = GetParam(argc, argv, &odometer_flag);
  if (param == (Param_t *)NULL) EXIT_ERROR("getting runtime parameters",
    "main");

  printf ("\nRunning lndcal ...\n");

  /* Validate the input metadata file */
  if (validate_xml_file (param->input_xml_file_name) != SUCCESS)
  {  /* Error messages already written */
      EXIT_ERROR("Failure validating XML file", "main");
  }

  /* Initialize the metadata structure */
  init_metadata_struct (&xml_metadata);

  /* Parse the metadata file into our internal metadata structure; also
     allocates space as needed for various pointers in the global and band
     metadata */
  if (parse_metadata (param->input_xml_file_name, &xml_metadata) != SUCCESS)
  {  /* Error messages already written */
    EXIT_ERROR("parsing XML file", "main");
  }

  /* Check to see if the gain and bias values were specified */
  if (!existRadGB (&xml_metadata))
    EXIT_ERROR("Gains and biases don't exist in XML file (TOA radiance gain "
      "and bias fields) for each band.  Make sure to utilize the latest LPGS "
      "MTL file for conversion to the ESPA internal raw binary format as the "
      "gains and biases should be in that file.", "main");
  
  /* Open input files */
  input = OpenInput (&xml_metadata);
  if (input == (Input_t *)NULL)
    EXIT_ERROR("setting up input from XML structure", "main");

  /* Get Lookup table */
  lut = GetLut(param, input->nband, input);
  if (lut == (Lut_t *)NULL) EXIT_ERROR("bad lut file", "main");

  nps6=  input->size_th.s;
  nls6=  input->size_th.l;
  nps =  input->size.s;
  nls =  input->size.l;

  /* Open the output files.  Raw binary band files will be be opened. */
  output = OpenOutput(&xml_metadata, input, param, lut, false /*not thermal*/);
  if (output == NULL) EXIT_ERROR("opening output file", "main");

  /* Allocate memory for the input buffer, enough for all reflectance bands */
  input_psize = sizeof(unsigned char);
  line_in = calloc (nps * nband_refl, input_psize);
   if (line_in == NULL) 
     EXIT_ERROR("allocating input line buffer", "main");

  /* Allocate memory for the input solar zenith representative band buffer */
  input_psize = sizeof(int16);
  line_in_sun_zen = calloc (nps, input_psize);
  if (line_in_sun_zen == NULL) 
    EXIT_ERROR("allocating input line buffer for solar zenith band", "main");

  /* Create and open output thermal band, if one exists */
  if ( input->nband_th > 0 ) {
    output_th = OpenOutput (&xml_metadata, input, param, lut, true /*thermal*/);
    if (output_th == NULL)
      EXIT_ERROR("opening output therm file", "main");

    /* Allocate memory for the thermal input and output buffer, only holds
       one band */
    line_out_th = calloc(nps6, sizeof(uint16_t));
    if (line_out_th == NULL) 
      EXIT_ERROR("allocating thermal output line buffer", "main");
  } else {
    printf("*** no output thermal file ***\n"); 
  }

  /* Allocate memory for a single output line for the image data */
  line_out = calloc (nps, sizeof (uint16_t));
  if (line_out == NULL) 
    EXIT_ERROR("allocating output line buffer", "main");

  /* Allocate memory for one line of the QA data. */
  line_qa = calloc (nps, sizeof(uint16_t));
  if (line_qa == NULL)
    EXIT_ERROR("allocating qa input line buffer", "main");

  /* Do for each THERMAL line */
  if (input->nband_th > 0) {
    for (iline = 0, curr_line = 0; iline < nls6; iline++, curr_line += nps6) {
      if (odometer_flag && (iline == 0 || iline == nls-1 || iline%100 == 0)) {
          printf("--- main loop BAND6 Line %d --- \r", iline);
          fflush(stdout);
      }

      /* Read the input thermal data */
      if (!GetInputLineTh(input, iline, line_in, &line_in[nps6]))
        EXIT_ERROR("reading input data for a line", "main");

      /* Read the input L1 QA data */
      if (!GetInputLineQA(input, iline, line_qa))
        EXIT_ERROR("reading input data for a line", "main");

      /* Handle the TOA brightness temp corrections */
      if (input-> nband_th == 1)
      {
        if (!Cal_TM_thermal(lut, input, line_in, line_out_th, line_qa, iline))
            EXIT_ERROR("doing calibration for a line", "main");
      }
      else if (input->nband_th == 2)
      {
        if (!Cal_ETM_thermal(lut, input, line_in, &line_in[nps6], line_out_th,
            line_qa, iline))
            EXIT_ERROR("doing calibration for a line", "main");
      }

      /* Write the results */
      ib=0;
      if (!PutOutputLine(output_th, ib, iline, line_out_th)) {
        sprintf(msgbuf,"write thermal error ib=%d iline=%d",ib,iline);
        EXIT_ERROR(msgbuf, "main");
      }
    } /* end loop for each thermal line */
  }
  if (odometer_flag) printf("\n");

  if (input->nband_th > 0)
    if (!CloseOutput(output_th))
      EXIT_ERROR("closing output thermal file", "main");

  /* Do for each REFLECTIVE line */
  for (iline = 0, curr_line = 0; iline < nls; iline++, curr_line+=nps) {
    if (odometer_flag && (iline == 0 || iline == nls-1 || iline%100 == 0)) {
        printf("--- main reflective loop Line %d ---\r", iline);
        fflush(stdout);
    }

    /* Read the input reflectance data */
    for (ib = 0; ib < input->nband; ib++) {
      if (!GetInputLine(input, ib, iline, &line_in[ib*nps]))
        EXIT_ERROR("reading input data for a line", "main");
    }
    
    /* Read input representative solar zenith band */
    if (!GetInputLineSunZen(input, iline, line_in_sun_zen))
      EXIT_ERROR("reading input solar zenith data for a line", "main");

    /* Read the input L1 QA data */
    if (!GetInputLineQA(input, iline, line_qa))
        EXIT_ERROR("reading input data for a line", "main");

    /* Handle the TOA reflectance corrections for every band, and write the
       results */
    for (ib = 0; ib < input->nband; ib++) {
      if (!Cal(param, lut, ib, input, &line_in[ib*nps], line_in_sun_zen,
        line_out, line_qa, iline))
        EXIT_ERROR("doing calibraton for a line", "main");

      if (!PutOutputLine(output, ib, iline, line_out))
        EXIT_ERROR("reading input data for a line", "main");
    } /* End loop for each band */
        
  } /* End loop for each line */

  if ( odometer_flag )printf("\n");

  /* Free the data arrays */
  free(line_out);
  line_out = NULL;
  free(line_in);
  line_in = NULL;
  free(line_in_sun_zen);
  line_in_sun_zen = NULL;
  free(line_out_th);
  line_out_th = NULL;
  free(line_qa);
  line_qa = NULL;

  /* Close input and output files */
  if (!CloseInput(input)) EXIT_ERROR("closing input file", "main");
  if (!CloseOutput(output)) EXIT_ERROR("closing input file", "main");

  /* Write the ENVI header for reflectance files */
  for (ib = 0; ib < output->nband; ib++) {
    /* Create the ENVI header file this band */
    if (create_envi_struct (&output->metadata.band[ib], &xml_metadata.global,
      &envi_hdr) != SUCCESS)
        EXIT_ERROR("Creating the ENVI header structure for this file.", "main");

    /* Write the ENVI header */
    strcpy (envi_file, output->metadata.band[ib].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        EXIT_ERROR("Writing the ENVI header file.", "main");
  }

  /* Write the ENVI header for thermal files */
  for (ib = 0; ib < output_th->nband; ib++) {
    /* Create the ENVI header file this band */
    if (create_envi_struct (&output_th->metadata.band[ib], &xml_metadata.global,
      &envi_hdr) != SUCCESS)
        EXIT_ERROR("Creating the ENVI header structure for this file.", "main");

    /* Write the ENVI header */
    strcpy (envi_file, output_th->metadata.band[ib].file_name);
    cptr = strchr (envi_file, '.');
    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
        EXIT_ERROR("Writing the ENVI header file.", "main");
  }

  /* Append the reflective and thermal bands to the XML file */
  if (append_metadata (output->nband, output->metadata.band,
    param->input_xml_file_name) != SUCCESS)
    EXIT_ERROR("appending reflectance and QA bands", "main");
  if (input->nband_th > 0) {
    if (append_metadata (output_th->nband, output_th->metadata.band,
      param->input_xml_file_name) != SUCCESS)
      EXIT_ERROR("appending thermal and QA bands", "main");
  }

  /* Free the metadata structure */
  free_metadata (&xml_metadata);

  /* Free memory */
  if (!FreeParam(param)) 
    EXIT_ERROR("freeing parameter stucture", "main");

  if (!FreeInput(input)) 
    EXIT_ERROR("freeing input file stucture", "main");

  if (!FreeLut(lut)) 
    EXIT_ERROR("freeing lut file stucture", "main");

  if (!FreeOutput(output)) 
    EXIT_ERROR("freeing output file stucture", "main");

  /* All done */
  printf ("lndcal complete.\n");
  return (EXIT_SUCCESS);
}
