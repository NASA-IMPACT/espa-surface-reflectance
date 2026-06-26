lndsr File Format Description
=============================

LEDAPS Product: lndsr
File Type: HDF
File Content Description: Landsat surface reflectance
Scientific Data Sets:
   band1 
   band2
   band3
   band4
   band5
   band7
   atmos_opacity
   fill_QA
   DDV_QA
   cloud_QA
   cloud_shadow_QA
   snow_QA
   land_water_QA
   adjacent_cloud_QA
   band6 


Essential SDSes
---------------
The essential SDS layers include surface reflectance for band 1-5 and 7 and quality flag. 
Each image data SDS layer is written as signed 16-bit integer type (short int) in 2-D data layer.  
 
SDS_name	long_name			valid_range		_FillValue	scale_factor 
band1		band 1 reflectance	-2000, 12000	-9999		0.0001
band2		band 2 reflectance	-2000, 12000    -9999		0.0001
band3		band 3 reflectance	-2000, 12000    -9999		0.0001
band4		band 4 reflectance	-2000, 12000    -9999		0.0001
band5		band 5 reflectance	-2000, 12000    -9999		0.0001
band7		band 7 reflectance	-2000, 12000    -9999		0.0001

Data Conversion
actual_reflectance = file_digital_number * scale_factor

Surface reflectance quality flags are saved as unsigned bytes in individual bands with values either on (255) or off (0). For example

SDS_name	long_name		valid_range
fill_QA		fill_QA			0, 255
cloud_QA	cloud_QA		0, 255


Optional SDSes 
--------------
The optional SDS layer include band 6 brightness temperature     
(after using lndapp or lndsrbm.ksh)
SDS_name     long_name		        unit	 valid_range	_FillValue	scale_factor 
band6	     band 6 temperature	    celsius	 -7000,7000		-9999		0.01

Data Conversion
actual_temperature = file_digital_number * scale_factor


Global Attributes
-----------------
Global attributes include metadata information parsed from input Landsat scene. An example
global metadata for a Landsat 7 ETM+ scene (p15r33) acquired on October 5, 2001 looks like:

        :DataProvider = "USGS/EROS" ;
        :Satellite = "LANDSAT_7" ;
        :Instrument = "ETM" ;
        :AcquisitionDate = "2001-10-05T17:42:00.000000Z" ;
        :Level1ProductionDate = "2002-12-30T00:00:00.000000Z" ;
        :SolarZenith = 48.210518f ;
        :SolarAzimuth = 152.53194f ;
        :WRS_System = "2" ;
        :WRS_Path = 15 ;
        :WRS_Row = 33 ;
        :NumberOfBands = '\6' ;
        :BandNumbers = '\1', '\2', '\3', '\4', '\5', '\7' ;
        :ShortName = "L7ESR" ;
        :LocalGranuleID = "L7ESR.a2001278.w2p015r033.020.2006356233850.hdf" ;
        :ProductionDate = "2006-12-22T23:38:50Z" ;
        :PGEVersion = "2.0.2" ;
        :ProcessVersion = "2.0.2" ;
        :WestBoundingCoordinate = -78.4317587424775 ;
        :EastBoundingCoordinate = -75.41061199149598 ;
        :NorthBoundingCoordinate = 39.94070883207912 ;
        :SouthBoundingCoordinate = 37.84141476946945 ;
        :ozone_source = "EARTHPROBE" ;
        :OrientationAngle = 0. ;
        :PixelSize = 28.5 ;
        :HDFEOSVersion = "HDFEOS_V2.4" ;
        :StructMetadata_0 = "GROUP=SwathStructure\n",
    "END_GROUP=SwathStructure\n",
    "GROUP=GridStructure\n",
    "\tGROUP=GRID_1\n",
    "\t\tGridName=\"Grid\"\n",
    "\t\tXDim=9029\n",
    "\t\tYDim=7993\n",
    "\t\tUpperLeftPointMtrs=(206568.000000,4421262.000000)\n",
    "\t\tLowerRightMtrs=(463894.500000,4193461.500000)\n",
    "\t\tProjection=GCTP_UTM\n",
    "\t\tZoneCode=18\n",
    "\t\tSphereCode=8\n",
    "\t\tGridOrigin=HDFE_GD_UL\n",
    "\t\tGROUP=Dimension\n",
    "\t\tEND_GROUP=Dimension\n",
    "\t\tGROUP=DataField\n",
    "\t\t\tOBJECT=DataField_0\n",
    "\t\t\t\tDataFieldName=\"band1\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_0\n",
    "\t\t\tOBJECT=DataField_1\n",
    "\t\t\t\tDataFieldName=\"band2\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_1\n",
    "\t\t\tOBJECT=DataField_2\n",
    "\t\t\t\tDataFieldName=\"band3\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_2\n",
    "\t\t\tOBJECT=DataField_3\n",
    "\t\t\t\tDataFieldName=\"band4\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_3\n",
    "\t\t\tOBJECT=DataField_4\n",
    "\t\t\t\tDataFieldName=\"band5\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_4\n",
    "\t\t\tOBJECT=DataField_5\n",
    "\t\t\t\tDataFieldName=\"band7\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_5\n",
    "\t\t\tOBJECT=DataField_6\n",
    "\t\t\t\tDataFieldName=\"aot_band1\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_6\n",
    "\t\t\tOBJECT=DataField_7\n",
    "\t\t\t\tDataFieldName=\"lndsr_QA\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_7\n",
    "\t\t\tOBJECT=DataField_8\n",
    "\t\t\t\tDataFieldName=\"nb_dark_pixels\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_8\n",
    "\t\t\tOBJECT=DataField_9\n",
    "\t\t\t\tDataFieldName=\"avg_dark_sr_b7\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_9\n",
    "\t\t\tOBJECT=DataField_10\n",
    "\t\t\t\tDataFieldName=\"std_dark_sr_b7\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_10\n",
    "\t\t\tOBJECT=DataField_12\n",
    "\t\t\t\tDataFieldName=\"band6\"\n",
    "\t\t\t\tDataType=DFNT_INT16\n",
    "\t\t\t\tDimList=(\"YDim\",\"XDim\")\n",
    "\t\t\tEND_OBJECT=DataField_12\n",
    "\t\tEND_GROUP=DataField\n",
    "\t\tGROUP=MergedFields\n",
    "\t\tEND_GROUP=MergedFields\n",
    "\tEND_GROUP=GRID_1\n",
    "END_GROUP=GridStructure\n",
    "GROUP=PointStructure\n",
    "END_GROUP=PointStructure\n",
    "END\n",
    "" ;
If lndsr data run through lndsrbm.ksh, additional line will be added 
    :Cloud_Mask_Algo_Version = "CMReflectanceBasedv1.0" ;


