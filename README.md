# 3DUS-US/MRI Rigid Registration
[Paper](https://link.springer.com/article/10.1007/s11548-023-02915-0)

We proposed a vessel-based 3DUS-CT/MRI rigid registration for liver tumor ablation. This [work](https://ieeexplore.ieee.org/abstract/document/9800921) shows the way of collecting 3D US images. To achieve the 3D US-CT/MRI registration task, firstly, we trained a nnUNet model to automatically segment the vessels from 3D US and CT/MRI data. Next, the vessel surface model (point clouds) and centerlines are extracted by using the [3D Slicer](https://www.slicer.org/) software. Lastly, the coarse-to-fine registration algorithm can help us align them together. The whole work is shown below.

<p align="center"><img src="1_updated.png" width="700" height="400"> </p>

## Citation
If you use this code for your research, please cite our paper:
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
```
## Introduction
