# -*- coding: utf-8 -*-
"""
Class for trajectory analysis of Amber MD simulation data
"""

import os
import sys

import pickle as pkl
import pytraj as pt
import pandas as pd
import numpy as np
import numpy.linalg as ln


class Protein:
    '''
    Class for trajectory analysis of Amber MD simulation data.
    Uses Pytraj to calculate some quantities of my choice.

    Example:
    >>> prot = cppy.Protein("1xwn_sh")
    
    Processing 1XWN
    
    PROTEIN INFORMATION
    PDB code: 1XWN
    residues: 166
    frames: <FRAMES>
    temperature: 300
    >>> print(prot.PDB)
    1XWN
    >>> # if data not found
    >>> prot = cppy.Protein("1xwn_sh", load_saved=True)
    
    Processing 1XWN
    Couldn't load data, processing new.
    
    PROTEIN INFORMATION
    PDB code: 1XWN
    residues: 166
    frames: <FRAMES>
    temperature: 300
    >>> # otherwise
    >>> prot = cppy.Protein("1xwn_sh", load_saved=True)
    
    Processing 1XWN
    Loaded data_1XWN
    
    PROTEIN INFORMATION
    PDB code: 1XWN
    residues: 166
    frames: <FRAMES>
    temperature: 300
    
    >>> # print radius of gyration
    >>> prot.rg
    15.068800201513314
    '''

    def __init__(self, prot, strd=0, mask="@CA", calc_all=False,
                 load_saved=False, save_new=False, temperature=300):
        self.prot = prot
        self.strd = strd
        self.mask = mask
        self.temperature = temperature
        # Protein code as it appears in RCSB Protein Data Bank
        self.PDB = prot[:-3].upper() if '_sh' in prot else prot.upper()
        self.data_path = 'temp{}K/data_{}'.format(self.temperature, self.PDB)
        self.absol = '..' #path to raw_data folder

        print('\nProcessing {}'.format(self.PDB))
        self.fetch_data(load_saved, calc_all, save_new)

    def fetch_data(self, load_saved, calc_all, save_new):
        '''Load saved data or trajectory file.
        If saved_data is True but no data is found it will load new data.
        '''
        if load_saved:
            try:
                self.load_data()
            except:
                # print(sys.exc_info())
                print("Couldn't load data, processing new.")
                self.prepare_all(calc_all, save_new)
        else:
            self.prepare_all(calc_all, save_new)

    def prepare_all(self, calc_all, save_new):
        self.prepare()
        if calc_all:
            self.calc_all()
        if save_new:
            self.save_data()

    def load_data(self):
        with open(os.path.join(self.absol, 'processed_data', self.data_path), 'rb') as f:
            # there's some issue with encoding and pickle in python 3
            loadedData = pkl.load(f, encoding='latin1')

            self.ppal = loadedData['ppal']
            self.rg = loadedData['radgyr']
            self.atoms = loadedData['resnum']
            self.dist = loadedData['distance']
            self.frames = loadedData['frames']
            self.cov = loadedData['covariance']
            self.xyz = loadedData['trajectory']
            self.corr = loadedData['correlation']
            self.SF_fluct = loadedData['SF_fluct']
            self.avg_dist = loadedData['avg_dist']
            self.corr_len = loadedData['correl_len']
            self.temperature = loadedData['temperature']
            self.fluct_norms = loadedData['fluct_norms']
            self.shape_factor = loadedData['shape_factor']
            self.fluctuations = loadedData['fluctuations']
            self.suceptibility = loadedData['suceptibility']
            self.avgd_corr = pd.DataFrame(loadedData['avgd_correlation'])

            print('Loaded data_{}'.format(self.PDB))
            print('\nPROTEIN INFORMATION')
            print('PDB code: {}'.format(self.PDB))
            print('residues: {}'.format(self.atoms))
            print('frames: {}'.format(self.frames))
            print('temperature: {}'.format(self.temperature))

    def prepare(self):
        '''Load trajectory, print information and
        perform initial alignment and fit.'''
        self.load_traj()
        self.atoms = self.traj.n_atoms
        self.frames = self.traj.n_frames

        print('\nPROTEIN INFORMATION')
        print('PDB code: {}'.format(self.PDB))
        print('residues: {}'.format(self.atoms))
        print('frames: {}'.format(self.frames))
        print('temperature: {}'.format(self.temperature))

        # Align molecule to principal axes
        pt.principal_axes(self.traj, dorotation=True)
        # Center molecule at the origin
        pt.center(self.traj, center='origin')
        # First rmsd fit
        pt.rmsd(self.traj)
        # Calculate average structure
        self.avg = pt.mean_structure(self.traj)
        # Perform rms fit to average structure
        pt.rmsd(self.traj, ref=self.avg)
        # Initialize arrays
        self.initialize_arrays()

    def initialize_arrays(self):
        # fluctuations
        self.fluct = []
        # magnitude of fluctuations
        self.fluct_norms = []
        # covariance coeficients matrix
        self.cov = np.zeros((self.traj.n_atoms, self.traj.n_atoms))
        # correlation coefficients matrix
        self.corr = np.zeros((self.traj.n_atoms, self.traj.n_atoms))
        # distance between residues matrix
        self.dist = np.zeros((self.traj.n_atoms, self.traj.n_atoms))
        # averaged correlations
        self.avgd_corr = pd.DataFrame()
        # correlation length
        self.corr_len = 0
        # shape factor fluctuations
        self.SF_fluct = []
        # suceptibility
        self.suceptibility = 0
        # time averaged MSD
        self.TA_MSD = []
        # autocorrelation
        self.autocorr = []

    def calc_all(self):
        # calculate fluctuations
        self.calc_fluct()
        # calculate covariances
        self.calc_cov()
        # calculate correlations
        self.calc_corr()
        # calculate principal axes (length)
        self.calc_ppal_axes()
        # calculate radius of gyrations
        self.rg = np.mean(pt.radgyr(self.traj))
        # computes average distance between residues (~3.8A = .38nm)
        tem_dist = [self.dist[i][i+1] for i in range(self.traj.n_atoms - 1)]
        self.avg_dist = np.mean(tem_dist)
        # calculate shape factor fluctuations and mean value
        self.shape_factor_fluct()
        self.shape_factor = np.mean(self.SF_fluct)
        # compute suceptibility
        self.calc_suceptibility()
        # compute TA-MSD
        # self.calc_TA_MSD() # not yet
        # compute autocorrelation in terms of the lag
        # self.calc_autocorr() # not yet

    def save_data(self):
        '''Save all data to binary file.'''
        dataset = {'PDB': self.PDB,
                   'ppal': self.ppal,
                   'radgyr': self.rg,
                   'resnum': self.atoms,
                   'distance': self.dist,
                   'frames': self.frames,
                   'covariance': self.cov,
                   'correlation': self.corr,
                   'SF_fluct': self.SF_fluct,
                   'avg_dist': self.avg_dist,
                   'fluctuations': self.fluct,
                   'trajectory': self.traj.xyz,
                   'correl_len': self.corr_len,
                   'temperature': self.temperature,
                   'fluct_norms': self.fluct_norms,
                   'shape_factor': self.shape_factor,
                   'suceptibility': self.suceptibility,
                   'avgd_correlation': self.avgd_corr.to_dict()}

        with open(os.path.join(self.absol, 'processed_data', self.data_path), 'wb') as f:
            pkl.dump(dataset, f)
            print('Saved data_{} in folder {}'.format(self.PDB, self.data_path))

    def load_traj(self):
        temp_dir = 'temp{}k'.format(self.temperature)
        # path to files
        files = os.path.join(self.absol, "raw_data", self.prot, temp_dir, 'input')
        # trajectory path
        traj = os.path.join(files, 'nowat.{}.nc'.format(self.prot))
        # topology path
        top = os.path.join(files, 'nowat.{}.prmtop'.format(self.prot))
        # load trajectory file
        self.traj = pt.load(traj, top, mask=self.mask, stride=self.strd)

    def avg_corr(self):
        '''Performs rolling average of correlation and calculated
        correlation length from the resulting curve.'''
        # flatten distance and correlation arrays, create dataframe
        flat_corr = np.ndarray.flatten(np.array(self.corr))
        flat_dist = np.ndarray.flatten(np.array(self.dist))
        corr_dict = {'distance': flat_dist, 'correlation': flat_corr}
        self.avgd_corr = pd.DataFrame(corr_dict)

        # clean data and perform rolling average
        self.avgd_corr.sort_values(by=['distance'], inplace=True)
        self.avgd_corr = self.avgd_corr[self.avgd_corr.distance != 0]
        self.avgd_corr = self.avgd_corr.rolling(self.atoms).mean()
        self.avgd_corr.dropna(axis=0, inplace=True)

        # correlation length defined as the first negative or zero correlation value
        cl_idx = (self.avgd_corr.correlation <= 0).idxmax()
        self.corr_len = self.avgd_corr.distance[cl_idx]
        print(u'\u2713 averaged correlation')

    def calc_fluct(self):
        '''Calculate atom/residue fluctuation from the average position.'''
        # Calculate fluctuations
        self.fluct = [frm.xyz - self.avg for frm in self.traj]
        self.fluct_norms = [[ln.norm(r) for r in diff] for diff in self.fluct]
        print(u'\u2713 calculated fluctuations')

    def calc_cov(self):
        '''Calculates covariance matrix from fluctuations.'''
        # Compute covariance coefficient matrix
        for i in range(self.atoms):
            for j in range(self.atoms):
                for frm in self.fluct:
                    self.cov[i][j] += np.sum([frm[i][q] * frm[j][q] for q in range(3)])
                self.cov[i][j] /= self.frames
        print(u'\u2713 calculated covariance')

    def calc_corr(self):
        '''Calculates correlation and distance matrices.'''
        # Compute correlation coefficient matrix
        # Calculate mutual distances between residues
        for i in range(self.atoms):
            for j in range(self.atoms):
                self.corr[i][j] = self.cov[i][j] / np.sqrt(self.cov[i][i] * self.cov[j][j])
                self.dist[i][j] = pt.distance(
                    self.traj, ':{} :{}'.format(i+1, j+1), frame_indices=[0])[0]
        print(u'\u2713 calculated correlations')

    def calc_ppal_axes(self):
        '''Calculate average of principal axes length.'''
        # Principal axes length
        xyz = self.traj.xyz
        frame_range = range(self.frames)

        minx = np.mean([min([q[0] for q in xyz[i]]) for i in frame_range])
        maxx = np.mean([max([q[0] for q in xyz[i]]) for i in frame_range])

        miny = np.mean([min([q[1] for q in xyz[i]]) for i in frame_range])
        maxy = np.mean([max([q[1] for q in xyz[i]]) for i in frame_range])

        minz = np.mean([min([q[2] for q in xyz[i]]) for i in frame_range])
        maxz = np.mean([max([q[2] for q in xyz[i]]) for i in frame_range])

        axisX = abs(maxx - minx) / 2
        axisY = abs(maxy - miny) / 2
        axisZ = abs(maxz - minz) / 2

        self.ppal = [axisX, axisY, axisZ]
        print(u'\u2713 calculated principal axes length')

    def shape_factor_fluct(self):
        '''Calculate shape factor time series.'''
        xyz = self.traj.xyz
        for frm in xyz:
            Lx = abs(max([q[0] for q in frm]) - min([q[0] for q in frm])) / 2
            Ly = abs(max([q[1] for q in frm]) - min([q[1] for q in frm])) / 2
            Lz = abs(max([q[2] for q in frm]) - min([q[2] for q in frm])) / 2
            S = self.atoms * (self.avg_dist ** 3) / (Lx * Ly * Lz)
            self.SF_fluct.append(S)
        print(u'\u2713 calculated shape factor')

    def calc_suceptibility(self):
        '''Calculate suceptibility as the sum of correlation values
        up to the correlation length.'''
        self.avg_corr()
        # x = self.avgd_corr['distance'][self.avgd_corr['distance'] <= self.corr_len]
        y = self.avgd_corr['correlation'][self.avgd_corr['distance'] <= self.corr_len]
        self.suceptibility = (self.shape_factor / self.atoms) * sum(y)  # ?
        print(u'\u2713 calculated suceptibility')

    def calc_autocorr(self):
        '''Autocorrelation.'''
        half = int(self.atoms / 2)
        R = pd.Series(pt.distance(
            self.traj, ':1-{} :{}-{}'.format(half, half, self.atoms)))
        # [R.autocorr(lag=i) for i in range(self.frames-5)]
        self.autocorr = pt.acorr(R)

    def calc_TA_MSD(self):
        '''Time averaged mean squared deviation.'''
        half = int(self.atoms / 2)
        dist = pt.distance(self.traj, ':1-{} :{}-{}'.format(half, half, self.atoms))
        aux = []
        for lag in range(len(dist) - 5):
            dif = np.mean([(dist[i + lag] - dist[i]) ** 2 for i in range(len(dist)-lag)])
            aux.append(dif)
        self.TA_MSD = aux
