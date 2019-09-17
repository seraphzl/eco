#include <stdio.h>
#include "gdal.h"
#include "cpl_conv.h" /* for CPLMalloc() */
#include "cpl_string.h"
#include <mkl.h>
#include <omp.h>

#define YBLOCKSIZE 1024
#define XBLOCKSIZE 53836
#define YSIZE 44947
#define XSIZE 53836
#define ALIGNMENT 32


int main(int argc, char const *argv[])
{
	/*define some constant*/
	char const *c3Filename = "c3.tif";
	char const *k6Filename = "k6.tif";
	char const *ls2Filename = "ls2.tif";
	char const *rFilename = "r.tif";
	char const *scFilename = "sc.tif";
	int const XOFF_c3 = 0, YOFF_c3 = 0;
	int const XOFF_k6 = 8252, YOFF_k6 = 1333;
	int const XOFF_ls2 = 574, YOFF_ls2 = 509;
	int const XOFF_r = 574, YOFF_r = 509;
	int xSize_c3, ySize_c3, count_c3, xBlockSize_c3, yBlockSize_c3;
	int nthreads = 4;

	// test
	// int yBlockSize = atoi(argv[1]);
	// nthreads = atoi(argv[2]);
	// char const *scFilename = argv[3];

    GDALAllRegister();

	/*Opening the File*/
    GDALDatasetH hDst_c3, hDst_k6, hDst_ls2, hDst_r;
    hDst_c3 = GDALOpen(c3Filename, GA_ReadOnly);
    hDst_k6 = GDALOpen(k6Filename, GA_ReadOnly);
    hDst_ls2 = GDALOpen(ls2Filename, GA_ReadOnly);
    hDst_r = GDALOpen(rFilename, GA_ReadOnly);
    if(hDst_c3 == NULL)
    {
        printf("Open %s ERROR\n", c3Filename);
        exit(1);
    }
    if(hDst_k6 == NULL)
    {
        printf("Open %s ERROR\n", k6Filename);
        exit(1);
    }
    if(hDst_ls2 == NULL)
    {
        printf("Open %s ERROR\n", ls2Filename);
        exit(1);
    }
    if(hDst_r == NULL)
    {
        printf("Open %s ERROR\n", rFilename);
        exit(1);
    }
    // printf("%d %s\n", argc, argv[1]);

    /*Getting Dataset Information*/
	// adfGeoTransform[0] /* top left x */
	// adfGeoTransform[1] /* w-e pixel resolution */
	// adfGeoTransform[2] /* 0 */
	// adfGeoTransform[3] /* top left y */
	// adfGeoTransform[4] /* 0 */
	// adfGeoTransform[5] /* n-s pixel resolution (negative value) */
	// get c3 info
	GDALDriverH hDriver_c3 = GDALGetDatasetDriver(hDst_c3);
	double adfGeoTransform_c3[6];
	GDALGetGeoTransform(hDst_c3, adfGeoTransform_c3);
	char const *projection_c3;
	projection_c3 = GDALGetProjectionRef(hDst_c3);
	xSize_c3 = GDALGetRasterXSize(hDst_c3);
	ySize_c3 = GDALGetRasterYSize(hDst_c3);
	count_c3 = GDALGetRasterCount(hDst_c3);

	// printf("Driver: %s/%s\n", GDALGetDriverShortName(hDriver), GDALGetDriverLongName(hDriver));
	// printf("Size: %d x %d x %d\n", xSize, ySize, count);
	// if(projection != NULL)
	// printf("Projection: '%s'\n", projection);
	// if(status_gt == CE_None)
	// {
    // printf("Origin = (%.6f, %.6f)\n", adfGeoTransform[0], adfGeoTransform[3]);
    // printf("Pixel Size = (%.6f, %.6f)\n", adfGeoTransform[1], adfGeoTransform[5]);
	// }

	/*Fetching a Raster Band*/
	GDALRasterBandH hBand_c3 = GDALGetRasterBand(hDst_c3, 1);
	GDALRasterBandH hBand_k6 = GDALGetRasterBand(hDst_k6, 1);
	GDALRasterBandH hBand_ls2 = GDALGetRasterBand(hDst_ls2, 1);
	GDALRasterBandH hBand_r = GDALGetRasterBand(hDst_r, 1);
	// int bGotMin, bGotMax;
	// double adfMinMax[2];
	GDALGetBlockSize(hBand_c3, &xBlockSize_c3, &yBlockSize_c3);
	printf("Block = %d x %d, Type = %s, ColorInterp = %s\n",
		   xBlockSize_c3, yBlockSize_c3,
		   GDALGetDataTypeName(GDALGetRasterDataType(hBand_c3)),
		   GDALGetColorInterpretationName(GDALGetRasterColorInterpretation(hBand_c3)));
	// adfMinMax[0] = GDALGetRasterMinimum( hBand, &bGotMin );
	// adfMinMax[1] = GDALGetRasterMaximum( hBand, &bGotMax );
	// if( ! (bGotMin && bGotMax) )
	//     GDALComputeRasterMinMax( hBand, TRUE, adfMinMax );
	// printf( "Min = %.3f, Max = %.3f\n", adfMinMax[0], adfMinMax[1] );
	// if( GDALGetOverviewCount(hBand) > 0 )
	//     printf( "Band has %d overviews.\n", GDALGetOverviewCount(hBand) );
	// if( GDALGetRasterColorTable( hBand ) != NULL )
	//     printf( "Band has a color table with %d entries.\n",
	//             GDALGetColorEntryCount(GDALGetRasterColorTable( hBand )));

	/*Reading Raster Data*/
	// 'k6': [8252, 1333, 53836, 44947]
	// int nDSXOff = 8252; //argv[1];
	// int nDSYOff = 1333; //argv[2];
	// int nXCount = 53836; //argv[3];
	// int nYCount = 10000; //argv[4];
	// int nBXSize = nXCount;
	// int nBYSize = nYCount;
	// short int *pafScanblock;
	// // int nXSize = GDALGetRasterBandXSize( hBand );
	// CPLErr status_read, status_write;
	// pafScanblock = (short int *) CPLMalloc(sizeof(short int)*nBYSize*nBXSize);
	// // status_read = GDALRasterIO( hBand, GF_Read, nDSXOff, nDSYOff, nXCount, nYCount,
	//             pafScanblock, nBXSize, nBYSize, GDT_Int16, 0, 0 );
	// if (status_read != CE_None)
	// {
	// 	printf("ERROR code : %d\n", status_read);
	// 	exit(1);		
	// }

	int BUFF_SIZE = YBLOCKSIZE*XBLOCKSIZE;
	// float *pafScanblock_c3, *pafScanblock_k6, *pafScanblock_ls2, *pafScanblock_r;
	// pafScanblock_c3 = (float *)mkl_malloc(sizeof(float)*BUFF_SIZE, ALIGNMENT);
	// pafScanblock_k6 = (float *)mkl_malloc(sizeof(float)*BUFF_SIZE, ALIGNMENT);
	// pafScanblock_ls2 = (float *)mkl_malloc(sizeof(float)*BUFF_SIZE, ALIGNMENT);
	// pafScanblock_r = (float *)mkl_malloc(sizeof(float)*BUFF_SIZE, ALIGNMENT);
	// short int *pafScanblock_buff = (short int *) CPLMalloc(sizeof(short int)*BUFF_SIZE);

	/*Writing Raster Data*/
	// char **papszOptions = NULL;
	// // papszOptions = CSLSetNameValue( papszOptions, "TILED", "YES" );
	// // papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "DEFLATE");
	// papszOptions = CSLSetNameValue(papszOptions, "BIGTIFF", "YES");
	// GDALDatasetH hDst_sc = GDALCreate(hDriver_c3, scFilename, XSIZE, YSIZE, 1, GDT_Float32, papszOptions);
	// GDALRasterBandH hBand_sc = GDALGetRasterBand(hDst_sc, 1);
	// GDALSetGeoTransform(hDst_sc, adfGeoTransform_c3);
	// GDALSetProjection(hDst_sc, projection_c3);
	// // set nodata value
	// double nd_c3 = GDALGetRasterNoDataValue(hBand_c3, NULL);
	// GDALSetRasterNoDataValue(hBand_sc, nd_c3);

	CPLErr status_read, status_write;
	int nthreads_max = omp_get_max_threads();
	printf("max omp threads: %d\n", nthreads_max);
	omp_set_num_threads(nthreads);
	printf("set working threads: %d\n", nthreads);
	printf("set blockSize: %d x %d\n", XBLOCKSIZE, YBLOCKSIZE);
	printf("xSize x ySize: %d x %d\n", XSIZE, YSIZE);

	// float *pafScanblock_sc;
	// pafScanblock_sc = (float *)mkl_malloc(sizeof(float)*BUFF_SIZE, ALIGNMENT);
	float sb_sc[YBLOCKSIZE*XBLOCKSIZE] __attribute__((aligned(ALIGNMENT)));
	float sb_c3[YBLOCKSIZE*XBLOCKSIZE] __attribute__((aligned(ALIGNMENT)));
	float sb_k6[YBLOCKSIZE*XBLOCKSIZE] __attribute__((aligned(ALIGNMENT)));
	float sb_ls2[YBLOCKSIZE*XBLOCKSIZE] __attribute__((aligned(ALIGNMENT)));
	float sb_r[YBLOCKSIZE*XBLOCKSIZE] __attribute__((aligned(ALIGNMENT)));
	// float *psb_sc __attribute__((aligned(ALIGNMENT)));
	// float *psb_c3 __attribute__((aligned(ALIGNMENT)));
	// float *psb_k6 __attribute__((aligned(ALIGNMENT)));
	// float *psb_ls2 __attribute__((aligned(ALIGNMENT)));
	// float *psb_r __attribute__((aligned(ALIGNMENT)));
	// psb_sc = sb_sc;
	// psb_c3 = sb_c3;
	// psb_k6 = sb_k6;
	// psb_ls2 = sb_ls2;
	// psb_r = sb_r;

	int yBlockSize = YBLOCKSIZE, xBlockSize = XBLOCKSIZE;

	// define struct to store the index of scanblock
	struct ScanblockIndex
	{
		float *psb_sc;
		float *psb_c3;
		float *psb_r;
		float *psb_k6;
		float *psb_ls2;
	} __attribute__((aligned(ALIGNMENT)));
	struct ScanblockIndex sbInd[BUFF_SIZE];
	struct ScanblockIndex *psbInd __attribute__((aligned(ALIGNMENT)));
	// build scanblock index
	psbInd = sbInd;
	for (int i = 0; i < BUFF_SIZE; i++)
	{
		psbInd->psb_sc = sb_sc+i;
		psbInd->psb_c3 = sb_c3+i;
		psbInd->psb_r = sb_r+i;
		psbInd->psb_k6 = sb_k6+i;
		psbInd->psb_ls2 = sb_ls2+i;
		psbInd++;
	}

	double s_initial = dsecnd();
	for (int i, yOff = 0; yOff < YSIZE; yOff += YBLOCKSIZE)
	{
		if (yOff+YBLOCKSIZE > YSIZE)
		{
			yBlockSize = YSIZE - yOff;
			BUFF_SIZE = yBlockSize*xBlockSize;
		}

		#pragma omp parallel
		{
			#pragma omp sections private(status_read)
			{
				#pragma omp section
				{
					status_read = GDALRasterIO(hBand_c3, GF_Read, XOFF_c3, YOFF_c3+yOff, xBlockSize, yBlockSize, sb_c3, xBlockSize, yBlockSize, GDT_Float32, 0, 0);
					if (status_read != CE_None)
						printf("Read c3 ERROR code: %d\n", status_read);
				}
				#pragma omp section
				{
					status_read = GDALRasterIO(hBand_r, GF_Read, XOFF_r, YOFF_r+yOff, xBlockSize, yBlockSize, sb_r, xBlockSize, yBlockSize, GDT_Float32, 0, 0);
					if (status_read != CE_None)
						printf("Read r ERROR code: %d\n", status_read);
				}
				#pragma omp section
				{
					status_read = GDALRasterIO(hBand_k6, GF_Read, XOFF_k6, YOFF_k6+yOff, xBlockSize, yBlockSize, sb_k6, xBlockSize, yBlockSize, GDT_Float32, 0, 0);
					if (status_read != CE_None)
						printf("Read k6 ERROR code: %d\n", status_read);
				}
				#pragma omp section
				{
					status_read = GDALRasterIO(hBand_ls2, GF_Read, XOFF_ls2, YOFF_ls2+yOff, xBlockSize, yBlockSize, sb_ls2, xBlockSize, yBlockSize, GDT_Float32, 0, 0);
					if (status_read != CE_None)
						printf("Read ls2 ERROR code: %d\n", status_read);
				}
			}

			// compute sc
			#pragma omp for simd private(i) schedule(static)
			for (i = 0; i < BUFF_SIZE; i++)
			{
				psbInd->psb_sc = sb_sc+i;
				*(psb_sc+i) = 1000.0f - *(psb_c3+i);
				*(psb_sc+i) *= *(psb_r+i);
				*(psb_sc+i) *= *(psb_k6+i);
				*(psb_sc+i) *= *(psb_ls2+i)*1.0e-11f;
				// sb_sc[i] = 1000.0f - sb_c3[i];
				// sb_sc[i] *= sb_r[i];
				// sb_sc[i] *= sb_k6[i];
				// sb_sc[i] *= sb_ls2[i]*1.0e-11f;
			}
		}
		// toc_loop = clock();
		// printf("Loop time - computation : %f seconds\n", (double)(toc_loop-tic_loop)/CLOCKS_PER_SEC);

		// double adfGeoTransform[6] = { 444720, 30, 0, 3751320, 0, -30 };
		// OGRSpatialReferenceH hSRS;
		// char *pszSRS_WKT = NULL;

		// hSRS = OSRNewSpatialReference( NULL );
		// OSRSetUTM( hSRS, 11, TRUE );
		// OSRSetWellKnownGeogCS( hSRS, "NAD27" );
		// OSRExportToWkt( hSRS, &pszSRS_WKT );
		// OSRDestroySpatialReference( hSRS );
		// CPLFree( pszSRS_WKT );

		// status_write = GDALRasterIO(hBand_sc, GF_Write, 0, yOff, xBlockSize, yBlockSize, sb_sc, xBlockSize, yBlockSize, GDT_Float32, 0, 0);
		// printf("Loop time - write raster : %f seconds\n", (double)(clock()-toc_loop)/CLOCKS_PER_SEC);
	}
	// clock_t toc = clock();
	// printf("Elapsed time : %f seconds\n", (double)(toc-tic)/CLOCKS_PER_SEC);
	double s_elapsed = dsecnd() - s_initial;
    printf("Elapsed time: %f seconds\n", s_elapsed);

	/*Closing the Files*/
	// CPLFree(pafScanblock_buff);
	// mkl_free(pafScanblock_sc);
	// mkl_free(pafScanblock_c3);
	// mkl_free(pafScanblock_k6);
	// mkl_free(pafScanblock_ls2);
	// mkl_free(pafScanblock_r);
	// GDALClose(hDst_sc);
	GDALClose(hDst_c3);
    GDALClose(hDst_k6);
    GDALClose(hDst_ls2);
    GDALClose(hDst_r);
	return 0;
}
