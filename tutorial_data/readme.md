# spm12w Tutorial Data Files

This directory contains all the necessary files for running through the spm12w tutorial. This allows you to run nearly all aspects of spm12w through it paces, ensuring that your installation is working. 

## About the provided datasets

All three tutorial subjects have signed a release form, allowing their anatomical and function data to be included in this tutorial. In order to maintain the anonymity of our tutorial subjects, anatomical
MRIs are not provided in the original PAR/REC. Instead, we provide NIFTI files which have been defaced, maintaining the brain but preventing 3D reconstruction of the face. In order to reduce the filesize of the example dataset, it is necessary that everyone convert the raw data from either PAR/REC to NIFTI or 
simply from archived raw nifti files (in the arch/nifti directory). Both approaches are automated with spm12w_prepare. We provide both PAR/REC and nifti files in order to allow Philips users to test conversion from raw PAR/REC files. Non Philips users can instead use the nifti files. See the tutorial documentation in the spm12w manual for more information. 

## Bugs and installation issues

The tutorial dataset has been thoroughly tested on multiple installations across platforms. If you run into issues running the dataset, ensure that spm12w has been properly installed (along with its dependencies) and that everything is in your matlab path. Further issues or bugs can be posted [here](https://github.com/ddwagner/spm12w/issues/new) or you can email me at [Dylan Wagner](mailto:dylan.d.wagner@gmail.com)