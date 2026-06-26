#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "hdf.h"
#include "mfhdf.h"

int verbose;
int copy_sds(int32 sdin_id,char *sds_name,int32 first_time_index,int32 sdout_id);

int main(int argc,char **argv) {
	int32 sdin_id,sdout_id,sds_id,dim_id,hdf_ind;
	int32 nb_datasets,nb_globattrs;
	int32 status,index,count,dims[MAX_VAR_DIMS],start[MAX_VAR_DIMS],edge[MAX_VAR_DIMS];
	int32 attr_index;
	int32 rank,data_type,attributes;
	int32 dim_index,dim_sizes[MAX_VAR_DIMS];
	int32 first_time_index;
	int32 sds_r1,sds_r2,sds_r3,sds_t3,sds_t4,sds_qc;
	int32 sds_lat,sds_lon,sds_ts,sds_tv,sds_phi;
	float32 ul_lat,ul_lon,lr_lat,lr_lon;
	char name[MAX_NC_NAME];
	VOIDP buffer;
	char label[200], unit[200], format[200], coordsys[200];
	int i,write_metadata;
	float32 time[4]={0.,6.,12.,18.};
	int16 doy;

   int32  num_lones,vgroup_id;  
   int32 *ref_array;    
   char   vgroup_name[VGNAMELENMAX], vgroup_class[VGNAMELENMAX];	
	char list_of_sds[4][10]={"pres","pr_wtr"};

	if (argc!=4) {
		fprintf(stderr,"usage: %s <input> <output> <doy>\n",argv[0]);
		exit(-1);
	}
	verbose=0;
/****
	open input
****/
	if ((sdin_id=SDstart(argv[1], DFACC_RDONLY))<0) {
   	fprintf(stderr,"can't open file %s\n",argv[1]);
		exit(-1);
	}
/****
	check and open output
****/
	write_metadata=0;
	if ((sdout_id=SDstart(argv[2], DFACC_RDONLY))<0) {
		if ((sdout_id=SDstart(argv[2], DFACC_CREATE))<0) {
   		fprintf(stderr,"can't create output %s\n",argv[2]);
			exit(-1);
		}
		write_metadata=1;
	} else {
		SDend(sdout_id);
		if ((sdout_id=SDstart(argv[2], DFACC_WRITE))<0) {
   		fprintf(stderr,"can't open output %s\n",argv[2]);
			exit(-1);
		}
	}
	doy=(int16)atoi(argv[3]);		
	
/****
	Determine the contents of the file
****/
   status = SDfileinfo(sdin_id, &nb_datasets, &nb_globattrs);
	if (verbose)
		printf ("Nb of global attributes = %d\n",nb_globattrs);
 	if (write_metadata) {
/****
		Copy global attributes to output if creating a new file
****/
  		for (index = 0; index < nb_globattrs; index++) {
     			status = SDattrinfo(sdin_id, index, name, &data_type, &count);
			if (verbose)
				printf ("attribute #%2d = %s (%d) size = %d\n",index,name,data_type,count);
			/* Allocate a buffer to hold the attribute data. */
			buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
			/* Read the attribute data. */
			status = SDreadattr(sdin_id, index, buffer);
			/* Write the attribute data. */
			status = SDsetattr(sdout_id, name, data_type,count,buffer);
			/* Free buffer */
			free(buffer);
  		}
/****
		Write Day Of Year to output
****/
		status = SDsetattr(sdout_id, "Day Of Year", DFNT_INT16,DFKNTsize(DFNT_INT16),&doy);
		
/****
		Copy lat/lon SDS
****/
		if (copy_sds(sdin_id,"lat",-1,sdout_id)) {
			fprintf(stderr,"ERROR: couldn't copy SDS lat ... ABORT\n");
			exit(-1);
		}
		if (copy_sds(sdin_id,"lon",-1,sdout_id)) {
			fprintf(stderr,"ERROR: couldn't copy SDS lon ... ABORT\n");
			exit(-1);
		}		
	}


	if (verbose)  	
		printf ("Nb of SDS = %d\n",nb_datasets);
/****
	Loop over SDS, get name and nb of attributes for each one
****/
   	for (index = 0; index < nb_datasets; index++) {
		sds_id = SDselect(sdin_id, index);
		status = SDgetinfo(sds_id, name, &rank, dim_sizes, \
                         &data_type, &attributes);

		if ((!strcmp(name,"pres"))||(!strcmp(name,"pr_wtr"))) {
			if (verbose) {
      				printf("SDS #%2d = %s  Rank = %d Type=%d  Nb Attrs=%d\n", index,name,rank,data_type,attributes);
				printf("\tDims = ");
				for (i=0;i<rank;i++) {
					if (i==0)
						printf("%d",dim_sizes[i]);
					else
						printf("x%d",dim_sizes[i]);
				}
			}
			first_time_index=(doy-1)*4;
			if (copy_sds(sdin_id,name,first_time_index,sdout_id)) {
				fprintf(stderr,"ERROR: couldn't copy SDS %s ... ABORT\n",name);
				exit(-1);
			}		


		}
	}
/****
		Close input & output
****/
	status = SDend(sdin_id);
	status = SDend(sdout_id);

}   



int copy_sds(int32 sdin_id,char *sds_name,int32 first_time_index,int32 sdout_id) {
	int32 index,sds_id,sdsout_id,status,start[MAX_VAR_DIMS],edge[MAX_VAR_DIMS];
	char name[MAX_NC_NAME];
	int32 rank,data_type,attributes,count,dim_id,dimout_id;
	int32 dim_index,dim_sizes[MAX_VAR_DIMS];
	VOIDP buffer;
	int buf_size,i;

			
        if ((index=SDnametoindex(sdin_id,sds_name))<0)
		return -1;
	sds_id = SDselect(sdin_id, index);
	status = SDgetinfo(sds_id, name, &rank, dim_sizes,&data_type, &attributes);
	if ((first_time_index >= 0) && (rank != 3))
		return -1;
	if ( rank==3 ) dim_sizes[0]=4;
	if ((sdsout_id=SDcreate(sdout_id,name,data_type,rank,dim_sizes))<0) 
		return -1;
	
	if (first_time_index < 0) {
		buf_size=1;
		for (i=0;i<rank;i++) {
			buf_size *= dim_sizes[i];
			start[i]=0;
			edge[i]=dim_sizes[i];
		}
		buf_size *= DFKNTsize(data_type);
		if ((buffer=(VOIDP)malloc(buf_size))==NULL)
			return -1;
		if (SDreaddata(sds_id,start,NULL,edge,buffer)<0)
			return -1;
	} else {
		start[0]= first_time_index/4;
		edge[0]=4;
		start[1]=0;
		edge[1]=dim_sizes[1];
		start[2]=0;
		edge[2]=dim_sizes[2];
		buf_size=4*dim_sizes[1]*dim_sizes[2]*DFKNTsize(data_type);
		if ((buffer=(VOIDP)malloc(buf_size))==NULL)
			return -1;
		if (SDreaddata(sds_id,start,NULL,edge,buffer)<0)
			return -1;
	}
	start[0]=0;
	if (SDwritedata(sdsout_id,start,NULL,edge,buffer)<0)
		return -1;

/****
	Loop over SDS attributes and copy
****/
      	for (index=0;index<attributes;index++) {
		/* Get information about the data set attribute. */
		status = SDattrinfo(sds_id, index, name, &data_type, &count);
		/* Allocate a buffer to hold the attribute data. */
		buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
		/* Read attribute  */
		status = SDreadattr(sds_id, index, buffer);
		/* Write attribute */
		SDsetattr(sdsout_id, name, data_type, count, buffer);
		/* Free buffer */
		free(buffer);
	}
		
/****
	Loop over SDS dimensions
****/
	for (index=0;index<rank;index++) {
		dim_id = SDgetdimid(sds_id, index);
		status = SDdiminfo(dim_id, name, &count, &data_type, &attributes);
		dimout_id=SDgetdimid(sdsout_id, index);
		SDsetdimname(dimout_id, name); 
		if (verbose)
			printf ("    Dim #%2d = %s  size = %d	data type = %d  nb attrs = %d\n",dim_index,name,count,data_type,attributes);
		if ((data_type != 0)&&(count != 0)) {
			/* Allocate a buffer to hold scale data. */
			buffer = (VOIDP)malloc(count * DFKNTsize(data_type));
			/* Read scale data. */
			status = SDgetdimscale(dim_id, buffer);
			/* Write scale data. */
			SDsetdimscale(dimout_id, count, data_type, buffer); 
			/* Free buffer */
			free(buffer); 
		}
	}
/****
	Close SDS
****/
	status = SDendaccess(sds_id);
	status = SDendaccess(sdsout_id);
	return 0;
}

