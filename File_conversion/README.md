This folder contains codes for unification of all the tiff files takes using two systems (Leica ans ISS into one single multipage file)
In order to use this code you will need to add the functions_u folder to MATLAB path
Each data set should be copied in a folder and next all the folders be put in a single folder
please remember not to mix multi channel data in a single folder, in other words your data in each folder should be either multichannel or single channel
I have assumed that the data from ISS system are only single channel
This code is sensitive to naming For Leica system the name should be similar to 20210624_EITB-000US-H2BGFP_02_ch0_z0_ch0.tif with first ch representing channel and second ch representing data type (int, G or S)
Also for ISS system the data name should look like 20191112_13002ORG-H2BGFP_lowEGF_CRCmedia_DC_Ch1_Z15_Hr1.tif or 20191112_13002ORG-H2BGFP_lowEGF_CRCmedia_S_Ch1_Z10_Hr1.tif
Code is developed by Soheil Soltani, for further inquiries please email ssoltani@eitm.org
