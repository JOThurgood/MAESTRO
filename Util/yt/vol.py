#!/usr/bin/env python

import matplotlib
matplotlib.use('agg')

import sys

import yt
import numpy as np
from yt.visualization.volume_rendering.api import \
    Scene, \
    VolumeSource


# this is for the wdconvect problem

def doit(plotfile):

    ds = yt.load(plotfile)
    ds.periodicity = (True, True, True)

    field = ('boxlib', 'radial_velocity')
    ds._get_field_info(field).take_log = False
        
    sc = Scene()


    # add a volume: select a sphere
    center = (0, 0, 0)
    R = (5.e8, 'cm')

    dd = ds.sphere(center, R)

    vol = VolumeSource(dd, field=field)
    vol.use_ghost_zones = True

    sc.add_source(vol)


    # transfer function
    vals = [-5.e6, -2.5e6, -1.25e6, 1.25e6, 2.5e6, 5.e6]
    sigma = 3.e5

    tf =  yt.ColorTransferFunction((min(vals), max(vals)))

    tf.clear()
    cm = "coolwarm"
    for v in vals:
        tf.sample_colormap(v, sigma**2, colormap=cm) #, alpha=0.2)

    sc.get_source(0).transfer_function = tf

    cam = sc.add_camera(ds, lens_type="perspective")        
    cam.resolution = (1280, 720)
    cam.position = 1.5*ds.arr(np.array([5.e8, 5.e8, 5.e8]), 'cm')
    
    # look toward the center -- we are dealing with an octant
    center = ds.domain_left_edge
    normal = (center - cam.position)
    normal /= np.sqrt(normal.dot(normal))

    cam.switch_orientation(normal_vector=normal,
                           north_vector=[0., 0., 1.])
    cam.set_width(ds.domain_width)

    #sc.annotate_axes()
    #sc.annotate_domain(ds)

    sc.render()
    sc.save("subchandra_test.png", sigma_clip=6.0)
    sc.save_annotated("subchandra_test_annotated.png", 
                      text_annotate=[[(0.05, 0.05), 
                                      "t = {}".format(ds.current_time.d),
                                      dict(horizontalalignment="left")],
                                     [(0.5,0.95), 
                                      "Maestro simulation of He convection on a white dwarf",
                                      dict(color="y", fontsize="24",
                                           horizontalalignment="center")]])

if __name__ == "__main__":

    # Choose a field
    plotfile = ""


    try: plotfile = sys.argv[1]
    except: sys.exit("ERROR: no plotfile specified")

    doit(plotfile)


        
