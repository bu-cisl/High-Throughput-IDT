# High-Throughput-IDT
MATLAB implementation of High throughput intensity diffraction tomography

### Getting Started
This code is based on the paper "High throughput intensity diffraction tomography with a computational microscope" [PDF](https://www.osapublishing.org/boe/abstract.cfm?uri=boe-9-5-2130
). It implements closed-form inversion based on Tikhonov regularization to recover 3D phase and absorption from 2D intensity measurements.

### Running the code
MATLAB is required to run this code
Run IDT_example.m

### Dataset
Raw data, saved as 16-bit image files in **LED_algae_data** is required to use the provided code. Each data file is a cropped intensity image acquired by illuminating a single LED in the illuminaiton array.

### Abstract
We demonstrate a motion-free intensity diffraction tomography technique that enables direct inversion of 3D phase and absorption from intensity-only measurements for weakly scattering samples.   We derive a novel linear forward model, featuring slice-wise phase and absorption transfer functions using angled illumination.  This new framework facilitates flexible and efficient data acquisition, enabling arbitrary sampling of the illumination angles.   The reconstruction algorithm performs 3D synthetic aperture using a robust, computation and memory efficient  slice-wise deconvolution to achieve resolution up to the incoherent limit.    We demonstrate our technique with thick biological samples having both sparse 3D structures and dense cell clusters.  We further investigate the limitation of our technique when imaging strongly scattering samples.  Imaging performance and the influence of multiple scattering is evaluated using a 3D sample consisting of stacked phase and absorption resolution targets.   This computational microscopy system is directly built on a standard commercial microscope with a simple LED array source add-on, and promises broad applications by leveraging the ubiquitous microscopy platforms with minimal hardware modifications.

### System overview
Intensity diffraction tomography from angled illumination. (a) The setup consists of a standard microscope with an LED array that allows flexible patterning of illumination angles. (b) Images are taken by varying the illumination angle. Each intensity spectrum of the raw data exhibits two shifted circles, whose shift is set by the illumination angle. (c) Corresponding phase (imaginary part) and amplitude (real part) transfer functions (TF) for the same set of illumination angles are visualized at various sample depths. (d) The slice-wise deconvolution algorithm outputs two 3D stacks, corresponding to the phase and absorption reconstruction.

![Alt Text](/images/figure_1.png)

### Results
Reconstruction of unstained MCF-7 cancer cells. (a) The full FOV PhC image (40x, 0.65NA).  (b) Phase reconstructions on a few cell regions, demonstrating its versatility and robustness in reconstructing both thin and thick samples.  (c) Phase reconstruction of a dense cell cluster across multiple slices.  The comparison with the physically scanned PhC images demonstrates that our IDT technique provides similar lateral resolution and axial sectioning capability.  (d)  Phase reconstruction of the cell clusters using symmetric and pseudorandom illumination patterns.  Our IDT framework allows flexibly designing the illumination pattern and the number of LEDs used.  The reconstruction algorithm produces high quality phase recovery as the number of images used is reduced, and remains robust even when the number of images is much fewer than the number of unknowns. 

![Alt Text](/images/figure_3.png)

## License
This project is licensed under the GNU General Public License v3- see the [LICENSE.md](LICENSE.md) file for details
