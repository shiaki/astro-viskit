#!/usr/bin/python

'''
  voronoi_plot.py: plot data on voronoi mesh.
  YJ Qin, May 2015, Shanghai
'''

import numpy as np
import scipy.spatial as spt

import matplotlib.patches     as patc
import matplotlib.collections as colle

def _cross_pt(xy_a, xy_b, rlim):

  xa, ya = xy_a; xb, yb = xy_b
  x1, x2, y1, y2 = rlim
  xc, yc = (x1 + x2) / 2., (y1 + y2) / 2.

  # parallel to X or Y?
  if (xa - xb) == 0.:
    if (ya - y1) * (yb - y1) <= 0.: return [(xa + xb) / 2., y1]
    else:                           return [(xa + xb) / 2., y2]
  if (ya - yb) == 0.:
    if (xa - x1) * (xb - x1) <= 0.: return [x1, (ya + yb) / 2.]
    else:                           return [x2, (ya + yb) / 2.]

  # non-trivial cases
  
  # x-lim crossing
  if (xa - x1) * (xb - x1) <= 0:
    ptx = [x1, ya + (yb - ya) * (x1 - xa) / (xb - xa)]
  elif (xa - x2) * (xb - x2) <= 0:
    ptx = [x2, ya + (yb - ya) * (x2 - xa) / (xb - xa)]
  else: ptx = None

  # y-lim crossing
  if (ya - y1) * (yb - y1) <= 0:
    pty = [xa + (xb - xa) * (y1 - ya) / (yb - ya), y1]
  elif (ya - y2) * (yb - y2) <= 0:
    pty = [xa + (xb - xa) * (y2 - ya) / (yb - ya), y2]
  else: pty = None

  if   ptx != None and pty == None: return ptx
  elif ptx == None and pty != None: return pty
  elif ptx != None and pty != None:
    if ((ptx[0] - xc) ** 2 + (ptx[1] - yc) ** 2) > \
       ((pty[0] - xc) ** 2 + (pty[1] - yc) ** 2): return pty
    else: return ptx
  else: raise RuntimeError('Life sucks')


def voronoi_plot(X, Y,
                 value       = None,
                 color_map   = 'jet', # str or ScalarMappable
                 draw_ridges = None,
                 field_lim   = None,
                 vmin        = None,
                 vmax        = None):

  vor = spt.Voronoi(np.vstack((X, Y)).T)
  vd_pts = np.array(vor.points)
  vd_vts = np.array(vor.vertices)
  rg_pts = vor.ridge_points
  rg_vts = vor.ridge_vertices
  vd_ptr = np.array(vor.point_region)

  id_valpt = np.zeros(len(vor.regions))
  for i_datapt, i_region in enumerate(vd_ptr):
    id_valpt[i_region] = i_datapt

  # check insiders
  if field_lim == None: field_lim = [X.min(), X.max(), Y.min(), Y.max()]
  id_ins = (field_lim[0] < vd_vts[:, 0]) * (vd_vts[:, 0] < field_lim[1]) * \
           (field_lim[2] < vd_vts[:, 1]) * (vd_vts[:, 1] < field_lim[3])
  is_inside = lambda i: ((i != -1) and id_ins[i])

  # check params
  if value == None: value = np.ones(X.size)

  # polygons
  poly, vals = [], []

  # iterate over polygons
  for i_region, poly_i in enumerate(vor.regions):

    # skip the first empty one
    if len(poly_i) == 0: continue

    # ridge poitns are outside the box, pass
    poly_ivt = poly_i + [poly_i[0]]
    id_ins_i = np.asarray(map(is_inside, poly_i))
    if np.all(~id_ins_i): continue

    # ridge pts are all inside, trivial case
    elif np.all(id_ins_i):
      poly_ixy = [vd_vts[ipt].tolist() for ipt in poly_i]

    # else, border-crossing one
    else:

      # for this polygon, xy coordinates
      poly_ixy = []

      # iterate over segments
      for ip_a, ip_b in zip(poly_ivt[:-1], poly_ivt[1:]):

        # ridge inside the cell, add point to the polygon
        if is_inside(ip_a) and is_inside(ip_b):
          poly_ixy.append(vd_vts[ip_a].tolist())

        else:

          ctr = vd_pts.mean(axis = 0)
          rlm = vd_pts.ptp(axis = 0)

          # running out of the box, insert two points
          if is_inside(ip_a) and ~is_inside(ip_b):
            #print 'A', ip_a, ip_b,
            poly_ixy.append(vd_vts[ip_a].tolist())

            # then deal with the second point
            if ip_b != -1:
              poly_ixy.append(_cross_pt(vd_vts[ip_a], vd_vts[ip_b],
                                        field_lim))
              #print ''
            else:
              # find the position of the far point (ip_b)
              rg_idx = rg_vts.index([-1, ip_a])
              #print rg_idx, rg_pts[rg_idx]
              tg  =  vd_pts[rg_pts[rg_idx][1]] - vd_pts[rg_pts[rg_idx][0]]
              tg /=  np.linalg.norm(tg)
              nm  =  np.array([-tg[1], tg[0]])
              mp  = (vd_pts[rg_pts[rg_idx][1]] + \
                     vd_pts[rg_pts[rg_idx][0]]) / 2.
              dr  =  np.sign(np.dot(mp - ctr, nm)) * nm
              fp  =  vd_vts[ip_a] + dr * rlm.max() 
              poly_ixy.append(_cross_pt(vd_vts[ip_a], fp, field_lim))

          # running into the box, insert one point
          elif ~is_inside(ip_a) and is_inside(ip_b):
            #print 'B', ip_a, ip_b,
            if ip_a != -1:
              poly_ixy.append(_cross_pt(vd_vts[ip_a], vd_vts[ip_b],
                                      field_lim))
              #print ''
            else:
              rg_idx = rg_vts.index([-1, ip_b])
              #print rg_idx, rg_pts[rg_idx]
              tg  =  vd_pts[rg_pts[rg_idx][1]] - vd_pts[rg_pts[rg_idx][0]]
              tg /=  np.linalg.norm(tg)
              nm  =  np.array([-tg[1], tg[0]])
              mp  = (vd_pts[rg_pts[rg_idx][1]] + \
                     vd_pts[rg_pts[rg_idx][0]]) / 2.
              dr  =  np.sign(np.dot(mp - ctr, nm)) * nm
              fp  =  vd_vts[ip_b] + dr * rlm.max() 
              poly_ixy.append(_cross_pt(vd_vts[ip_b], fp, field_lim))

          else: # none of them are inside the box, pass
            pass

      # check if hetro-cross 
      poly_ixy_c, poly_ixy_t = poly_ixy + [poly_ixy[0]], []
      for ivt in xrange(len(poly_ixy)):
        poly_ixy_t.append(poly_ixy_c[ivt])
        # lower y, left x
        if (poly_ixy_c[ivt][0]     == field_lim[0] and \
            poly_ixy_c[ivt + 1][1] == field_lim[2]) or \
           (poly_ixy_c[ivt + 1][0] == field_lim[0] and \
            poly_ixy_c[ivt][1]     == field_lim[2]):
          poly_ixy_t.append([field_lim[0], field_lim[2]])
        # upper y, left x
        if (poly_ixy_c[ivt][0]     == field_lim[0] and \
            poly_ixy_c[ivt + 1][1] == field_lim[3]) or \
           (poly_ixy_c[ivt + 1][0] == field_lim[0] and \
            poly_ixy_c[ivt][1]     == field_lim[3]):
          poly_ixy_t.append([field_lim[0], field_lim[3]])
        # lower y, right x
        if (poly_ixy_c[ivt][0]     == field_lim[1] and \
            poly_ixy_c[ivt + 1][1] == field_lim[2]) or \
           (poly_ixy_c[ivt + 1][0] == field_lim[1] and \
            poly_ixy_c[ivt][1]     == field_lim[2]):
          poly_ixy_t.append([field_lim[1], field_lim[2]])
        # upper y, right x
        if (poly_ixy_c[ivt][0]     == field_lim[1] and \
            poly_ixy_c[ivt + 1][1] == field_lim[3]) or \
           (poly_ixy_c[ivt + 1][0] == field_lim[1] and \
            poly_ixy_c[ivt][1]     == field_lim[3]):
          poly_ixy_t.append([field_lim[1], field_lim[3]])
      poly_ixy = poly_ixy_t

    # add the region into a list for plotting
    if poly_ixy != []:
      poly.append(poly_ixy)
      vals.append(value[id_valpt[i_region]])
    
  # Now turn poly into patch collection
  if (not draw_ridges): ec_param = 'none'
  else: ec_param = None

  poly_patches = [patc.Polygon(poly_i, True, ec = ec_param) for poly_i in poly]
  poly_collec  = colle.PatchCollection(poly_patches, cmap = color_map, edgecolor = ec_param)
  poly_collec.set_array(np.array(vals))
  poly_collec.set_clim(vmin, vmax)

  return poly_collec

# test
if False:

  ins = lambda X: abs(X[0]) < 1. and abs(X[1]) < 1.

  N = 20
  x1 = 2.4 * (np.random.random(N) - 0.5)
  x2 = 2.4 * (np.random.random(N) - 0.5)
  y1 = 2.4 * (np.random.random(N) - 0.5)
  y2 = 2.4 * (np.random.random(N) - 0.5)

  p1 = np.vstack((x1, y1)).T
  p2 = np.vstack((x2, y2)).T

  p1i, p2i = [], []
  for i in range(N):
    if ins(p1[i]) != ins(p2[i]):
      p1i.append(p1[i])
      p2i.append(p2[i])

  print len(p1i), len(p2i)

  import matplotlib.pyplot as plt
  fig = plt.figure()
  ax  = fig.add_subplot(1, 1, 1)
  for i in range(len(p1i)):
    ax.plot([p1i[i][0], p2i[i][0]], [p1i[i][1], p2i[i][1]])
    itp = _cross_pt(p1i[i], p2i[i], [-1., 1., -1., 1.])
    ax.scatter(itp[0], itp[1])
  plt.plot([-1, -1], [-1, 1], c = 'k')
  plt.plot([1, 1], [-1, 1], c = 'k')
  plt.plot([-1, 1], [-1, -1], c = 'k')
  plt.plot([-1, 1], [1, 1], c = 'k')

  plt.show()

if False:

  x = np.random.randn(50)
  y = np.random.randn(50)

  cc = voronoi_plot(x, y, field_lim = [-2.5, 2.5, -2.5, 2.5])

  import matplotlib.pyplot as plt

  fig = plt.figure(figsize = (15., 6.))
  ax = fig.add_subplot(1, 2, 1)
  ax2 = fig.add_subplot(1, 2, 2)
  vd = spt.Voronoi(np.vstack((x, y)).T)
  spt.voronoi_plot_2d(vd, ax2)
  ax.set_xlim(-3., 3.)
  ax.set_ylim(-3., 3.)
  ax.axvline(x = -2.5, linestyle = 'dashed')
  ax.axvline(x = 2.5, linestyle = 'dashed')
  ax.axhline(y = -2.5, linestyle = 'dashed')
  ax.axhline(y = 2.5, linestyle = 'dashed')
  ax2.set_xlim(-3., 3.)
  ax2.set_ylim(-3., 3.)
  ax2.axvline(x = -2.5, linestyle = 'dashed')
  ax2.axvline(x = 2.5, linestyle = 'dashed')
  ax2.axhline(y = -2.5, linestyle = 'dashed')
  ax2.axhline(y = 2.5, linestyle = 'dashed')

  for icc in cc:

    ict = icc + [icc[0]]
    for i in range(len(icc)):
      ax.plot([ict[i][0], ict[i + 1][0]],
              [ict[i][1], ict[i + 1][1]], c = 'k')

  for ivt, vtx in enumerate(vd.vertices):
    ax.annotate('%u'%ivt, xy = vtx)

  plt.show()

if False:

  x = np.random.randn(120)
  y = np.random.randn(120)
  z = np.random.randn(120)

  import matplotlib.pyplot as plt

  fig = plt.figure(figsize = (8., 6.))
  ax = fig.add_subplot(1, 1, 1)

  pat = voronoi_plot(x, y, z, color_map = 'gray_r',field_lim = [-2., 2., -2., 2.])
  ax.set_xlim(-2., 2.)
  ax.set_ylim(-2., 2.)
  ax.add_collection(pat)

  plt.show()
