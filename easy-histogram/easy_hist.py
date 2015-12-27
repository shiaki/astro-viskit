#!/usr/bin/python

'''
  adapt_hist.py: Create 2d/3d histogram with function evaluation.
  YJ Qin, Apr. 2015 / Shanghai
'''

import numpy as np
import itertools as itt

class easy_hist2d:

  def __init__(self, X, Y,
               hist_bins = (32, 32),
               hist_lim  = 'auto'):

    self.N_pts, self.N_bins = 0, None
    self.X_ax, self.Y_ax, self.X_id, self.Y_id = None, None, None, None
    self.hist = []

    if hist_lim == 'auto': hist_lim = [[X.min(), X.max()], [Y.min(), Y.max()]]

    N_pts = X.size
    pt_id = np.arange(N_pts)

    x_ax = np.linspace(hist_lim[0][0], hist_lim[0][1], hist_bins[0] + 1)
    y_ax = np.linspace(hist_lim[1][0], hist_lim[1][1], hist_bins[1] + 1)

    cid_x = np.digitize(X, x_ax)
    cid_y = np.digitize(Y, y_ax)

    # put them into cells
    self.hist = []
    for ic_x, ic_y in itt.product(xrange(1, hist_bins[0] + 1),
                                  xrange(1, hist_bins[1] + 1)):
      if ic_y == 1: self.hist.append([])
      id_msk = (cid_x == ic_x) * (cid_y == ic_y)
      id_ins = np.extract(id_msk, pt_id)
      self.hist[-1].append(id_ins)

    # set attr
    self.N_pts  = N_pts
    self.N_bins = hist_bins

    self.X_ax, self.Y_ax = x_ax, y_ax
    self.X_id, self.Y_id = cid_x, cid_y

  '''
  def __del__(self):
    self.N_pts, self.N_bins = 0, None
    self.X_ax, self.Y_ax, self.X_id, self.Y_id = None, None, None, None
    self.hist = []
  '''

  def evaluate(self, fc, args):

    #print 'hist size:', len(self.hist)

    rval = []
    for ic_x, ic_y in itt.product(xrange(self.N_bins[0]),
                                  xrange(self.N_bins[1])):
      if ic_y == 0: rval.append([])
      rval[-1].append(fc(self.hist[ic_x][ic_y], args))
      #if ic_x == 0 and ic_y == 0: print self.hist[ic_x][ic_y][:5] # DBG

    # unpack them
    if isinstance(rval[0][0], (tuple, list)):
      rv_lst = []
      for i_rv in range(len(rval[0][0])):
        rv_lst.append([])
        for ic_x, ic_y in itt.product(xrange(self.N_bins[0]),
                                      xrange(self.N_bins[1])):
          if ic_y == 0: rv_lst[-1].append([])
          rv_lst[-1][-1].append(rval[ic_x][ic_y][i_rv])
      return tuple([np.array(A) for A in rv_lst])

    else: return np.array(rval)
