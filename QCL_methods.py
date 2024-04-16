import numpy as np
from numpy.lib import scimath
from numpy.fft import fftshift, ifftshift, fft2, ifft2, fftfreq
from math import pi, floor, ceil
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.transforms import Bbox
import matplotlib.patches as patches
import h5py
import hdf5storage
import time
import gc 
import os
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.abspath('')))

class QCLFigure:

    def __init__(self, nrows = 1, ncols = 1, width = 5, height = 4, wspace = 0.2, hspace = 0.2):  

      self.fig, self.axes = plt.subplots(nrows = nrows, ncols = ncols)
      self.fig.set_size_inches(width, height)
      self.fig.subplots_adjust(wspace = wspace, hspace = hspace)

      if nrows is 1 and ncols is 1:
         self.figureLayout = 'Single'
      elif nrows is 1 or ncols is 1:
         self.figureLayout = '1D'
      else:    
         self.figureLayout = '2D'

    def plot(self, field, pixelSize, pos = (0,0), title = '', zoom = 1, 
             plot_lens_aperture = False, lens_aperture_radius = 1, plot_center_circle = False, circle_radius = 1, plot_integration_area = False,
             colormap = 'hot', plot_colorbar = False, animated = False):

        if self.figureLayout is 'Single':
            self.plot_axis = self.axes
        elif self.figureLayout is '1D':
            self.plot_axis = self.axes[pos[0] + pos[1]]
        elif self.figureLayout is '2D': 
            self.plot_axis = self.axes[pos[0], pos[1]]


        dimX = field.shape[0] 
        field_crop = arrayCenterCrop(field, dimX//zoom)
        cdimX, cdimY = field_crop.shape
        
        self.plot_axis.clear()
        if plot_lens_aperture is True:
            lens_aperture_circle = plt.Circle((cdimX//2, cdimY//2), lens_aperture_radius//pixelSize, color = [1,1,1], fill = False)
            self.plot_axis.add_artist(lens_aperture_circle)

        if plot_center_circle is True:   
            center_spot = plt.Circle((cdimX//2, cdimY//2), circle_radius, color = [1,0,0], fill = False)
            self.plot_axis.add_artist(center_spot)

        if plot_integration_area is True:    
            QCL_rect = patches.Rectangle((cdimX//2 - (32 - 5), cdimY//2 - (32 - 5)), 63 - 2*5, 63 - 2*5, linewidth = 1, edgecolor= [0, 1, 0], facecolor='none')
            self.plot_axis.add_patch(QCL_rect) 

        self.axes_image = self.plot_axis.imshow(field_crop, cmap=plt.get_cmap(colormap), animated = animated)
        self.plot_axis.set_title(title)
        self.plot_axis.set(xlabel='x [mm]', ylabel='y [mm]')
        self.plot_axis.set_xticks(np.linspace(0, cdimX-1, 3))
        self.plot_axis.set_xticklabels([f'{-pixelSize*cdimX/2*1e3:.2f}', 0, f'{pixelSize*cdimX/2*1e3:.2f}'], color='black')
        self.plot_axis.set_yticks(np.linspace(0, cdimY-1, 3))
        self.plot_axis.set_yticklabels([f'{-pixelSize*cdimY/2*1e3:.2f}', 0, f'{pixelSize*cdimY/2*1e3:.2f}'], color='black')
        if plot_colorbar == True:
            self.fig.colorbar(self.axes_image, ax=self.plot_axis)        


def arrayCenterCrop(array,Ncrop):
    y,x = array.shape
    startx = x//2 - Ncrop//2
    starty = y//2 - Ncrop//2    
    return array[starty:starty + Ncrop, startx:startx + Ncrop] 

def propagate(E_input, propagator, distance):

    N = len(E_input)
    Sk = (fftshift(fft2(ifftshift(E_input)))/N**2).astype(np.complex64)
    E_output = (fftshift(ifft2(ifftshift(Sk*np.exp(1j*propagator*distance))))*N**2).astype(np.complex64)
    return E_output 

def back_propagate(E_input, propagator, distance):

    N = len(E_input)
    Sk = (fftshift(fft2(ifftshift(E_input)))/N**2).astype(np.complex64)
    E_output = (fftshift(ifft2(ifftshift(Sk*np.exp(-1j*propagator.conj()*distance))))*N**2).astype(np.complex64)
    return E_output

def mirror(E_input):
    N = len(E_input)
    Sk = (fftshift(fft2(ifftshift(E_input)))/N**2).astype(np.complex64)
    Sk = np.flipud(Sk)
    Sk = np.fliplr(Sk)
    E_output = (fftshift(ifft2(ifftshift(Sk)))*N**2).astype(np.complex64)
    return E_output

def mirror_phase(E_input):
    E_output = np.abs(E_input)*np.exp(1j*(np.angle(E_input) + np.pi))
    return E_output    

def applyQuadraticPhase(E_input, lens_focal_length, wavelength, R_squared):
    return (E_input*np.exp(-1j*(pi/(wavelength*lens_focal_length))*R_squared)).astype(np.complex64) 

def there_and_back_and_integrate(E_QCL_sim_input, QCL_offset_x, QCL_offset_y, sample_offset_Z, d1, d2, d3, f1, f2, propagator, lens_aperture, wavelength, R_squared, X, Y, dimX, dimY,
                                 knife_edge = False, edge_type = 'x', knife_pix_offset = 0):
   
    tic = time.time()
    # propagate QCL output by d1
    E_enter_L1_forward = propagate(E_QCL_sim_input, propagator, d1)

    # apply lens clear aperture
    E_enter_L1_forward_apertured = (E_enter_L1_forward*lens_aperture).astype(np.complex64)

    # apply quadratic phase of the lens
    E_exit_L1_forward = applyQuadraticPhase(E_enter_L1_forward_apertured, wavelength, f1, R_squared)



    # propagate  from lens L1 to L2 (distance d2)
    E_enter_L2_forward = propagate(E_exit_L1_forward, propagator, d2)

    # apply lens clear aperture
    E_enter_L2_forward_apertured = (E_enter_L2_forward*lens_aperture).astype(np.complex64)

    # apply quadratic phase of the lens     
    E_exit_L2_forward = applyQuadraticPhase(E_enter_L2_forward_apertured, wavelength, f2, R_squared)

    # propagate  from lens L2 to sample (distance d3 + sample_offset_Z)
    E_sample_forward = propagate(E_exit_L2_forward, propagator, d3 + sample_offset_Z)



    # sample is a mirror
    E_sample_backward = mirror_phase(E_sample_forward)

    # apply knife edge
    if knife_edge and edge_type is 'x':
        knife_edge = (X >= X[0,dimX//2 + knife_pix_offset]) # knife edge from -rangeX to (knife_pix_x*pixelSizeX)
        E_sample_backward = E_sample_backward*knife_edge 
    elif knife_edge and edge_type is 'y':    
        knife_edge = (Y >= Y[dimY//2 + knife_pix_offset, 0]) # knife edge from -rangeY to (knife_pix_y*pixelSizeY) 
        E_sample_backward = E_sample_backward*knife_edge 

    # back-propagate from sample to L2 (distance d3 + sample_offset_Z)
    E_enter_L2_backward = propagate(E_sample_backward, propagator, d3 + sample_offset_Z)

    # apply lens clear aperture
    E_enter_L2_backward_apertured = (E_enter_L2_backward*lens_aperture).astype(np.complex64)

    # apply quadratic phase of the lens
    E_exit_L2_backward = applyQuadraticPhase(E_enter_L2_backward_apertured, wavelength, f2, R_squared)



    # back-propagate from L2 to L1 (distance d2)
    E_enter_L1_backward = propagate(E_exit_L2_backward, propagator, d2)

    # apply lens clear aperture
    E_enter_L1_backward_apertured = (E_enter_L1_backward*lens_aperture).astype(np.complex64)

    # apply quadratic phase of the lens
    E_exit_L1_backward = applyQuadraticPhase(E_enter_L1_backward_apertured, wavelength, f1, R_squared)

    # back-propagate from L1 to QCL (distance d1)
    E_QCL_coupled = propagate(E_exit_L1_backward, propagator, d1)

    toc = time.time()
    print(f'Iteration time {toc-tic} seconds.')

    # Test git
    
    return E_QCL_coupled


                     

if __name__ is '__main__':
    print('The module can run')        