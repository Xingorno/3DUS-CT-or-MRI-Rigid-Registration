
<div align="center">
<h1>3DUS-US/MRI Rigid Registration</h1>
<a href="https://xingorno.github.io/3DLIVUS/"><img src="https://img.shields.io/badge/Project-red" alt="Project Page"></a>
<a href="https://link.springer.com/article/10.1007/s11548-023-02915-0" target="_blank" rel="noopener noreferrer">
  <img src="https://img.shields.io/badge/Publication-IPCAI2023-green" alt="Publication">
</a>
</div>

<p align="justify"> We proposed a vessel-based 3DUS-CT/MRI rigid registration for liver tumour ablation. Our <a href="https://ieeexplore.ieee.org/abstract/document/9800921"> previous work </a> shows the way of collecting 3D US images. To achieve the 3D US-CT/MRI registration task, firstly, we trained a nnUNet model to automatically segment the vessels from 3D US and CT/MRI data. Next, the vessel surface models (point clouds) and centerlines are extracted by using the <a href="https://www.slicer.org/"> 3D Slicer </a> software. Lastly, the coarse-to-fine registration algorithm can help us align them together. The worklow is shown below.
</p>

<p align="center"><img src="1_updated.png" width="700" height="400"> </p>

## Citation
If you use this code for your research, please cite our publications:
```
@article{xing20233d,
  title={3D US-CT/MRI registration for percutaneous focal liver tumor ablations},
  author={Xing, Shuwei and Romero, Joeana Cambranis and Roy, Priyanka and Cool, Derek W and Tessier, David and Chen, Elvis CS and Peters, Terry M and Fenster, Aaron},
  journal={International Journal of Computer Assisted Radiology and Surgery},
  volume={18},
  number={7},
  pages={1159--1166},
  year={2023},
  publisher={Springer}
}
@article{xing20223d,
  title={3d us-based evaluation and optimization of tumor coverage for us-guided percutaneous liver thermal ablation},
  author={Xing, Shuwei and Romero, Joeana Cambranis and Cool, Derek W and Mujoomdar, Amol and Chen, Elvis CS and Peters, Terry M and Fenster, Aaron},
  journal={IEEE Transactions on Medical Imaging},
  volume={41},
  number={11},
  pages={3344--3356},
  year={2022},
  publisher={IEEE}
}
```
## Introduction
- ./vessel_segmentation_nnUNet. The well-trained nnUNet model for segmenting vessels directly from 3D US images. The checkpoint can be directly used, combined with the configured [nnUNet](https://github.com/MIC-DKFZ/nnUNet) enviroment.
- main.m. This file includes the complete coarse-to-refine registration algorithm to align the 3D US and CT/MRI images. Note that the vessels should be formatted as surface models represented as surface point clouds, and centerlines represented as a bunch of sampled points. [CloudCompare](https://www.danielgm.net/cc/) can be used to extract the surface point clouds, and 3D Slicer is used for extracting [centerlines](https://github.com/vmtk/SlicerExtension-VMTK). Both exported file formats are .txts.
- Registration_evaluation_TRE/centerlineDistance.m. These two files are used for evaluating the registration accuracy. "ReadRegisteredCenterlineParameters.m" shows you how to read the registered parameters.

***(If you can not download the trained models due to the huge size, please go to file "/vessel_segmentation_nnUNet/dropbox_link.txt". The shared dropbox link is included, you can download that externally.)***

## Setup
- Matlab 2023b. We used Matlab 2023b to run the registration algorithm. For Matlab, it is quite flexible, any different versions can also be used to run this code.
- nnUNet environment. Please follow this [repo](https://github.com/MIC-DKFZ/nnUNet) to configure the environment, which is to use our trained vessel segmentation model.
- 3D slicer (VMTK module). This [module](https://github.com/vmtk/SlicerExtension-VMTK) is to automatically extract the vessel centerlines and export the compatible format for our registration algorithm.
- CloudCompare. This software is to convert the vessel surface models (.stl or .obj) to the point clouds format (.txt).
 
