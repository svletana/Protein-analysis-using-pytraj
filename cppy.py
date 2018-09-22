# -*- coding: utf-8 -*-"""
"""
Class for processing Amber MD simulation data
"""

from __future__ import division
import os
import pickle as pkl
import pytraj as pt
import pandas as pd
import numpy as np
import numpy.linalg as ln

pdbs300 = pd.read_pickle("pdbCodes300K.dat")  # temperature: 300K
pdbs350 = pd.read_pickle("pdbCodes350K.dat")  # temperature: 350K


class Protein:
    def __init__(self, prot, strd=0, process_new=False, save_new=False, temperature=300):

        self.prot = prot
        self.strd = strd
        self.temperature = temperature
        if '_sh' in prot:
            self.PDBcode = prot[:-3].upper()
        else:
            self.PDBcode = prot.upper()

        print('\nProcessing {}'.format(self.PDBcode))
        self.data_path = 'processedData/temp{}K/data_{}'.format(self.temperature, self.PDBcode)

        if not process_new:
            try:
                f = open(self.data_path, 'rb')
                loadedData = pkl.load(f)

                self.atoms = loadedData['resnum']
                self.frames = loadedData['frames']
                self.temperature = loadedData['temperature']
                self.rg = loadedData['radgyr']
                self.avg_dist = loadedData['avg_dist']
                self.shape_factor = loadedData['shape_factor']
                self.SF_fluct = loadedData['SF_fluct']
                self.suceptibility = loadedData['suceptibility']
                self.corr_len = loadedData['correl_len']
                self.fluct_norms = loadedData['fluct_norms']
                self.avgd_corr = pd.DataFrame(loadedData['avgd_correlation'])
                self.fluctuations = loadedData['fluctuations']
                self.cov = loadedData['covariance']
                self.corr = loadedData['correlation']
                self.dist = loadedData['distance']
                self.traj = loadedData['trajectory']
                self.ppal = loadedData['ppal']

                f.close()
                print('Loaded data_{}'.format(self.PDBcode))

            except IOError:
                print('File not found, processing new trajectory')
                process_new = True

        if process_new:
            self.loadtraj()
            self.atoms = self.traj.n_atoms
            self.frames = self.traj.n_frames
            print('\nPROTEIN INFORMATION')
            print('PDB code: {}'.format(self.PDBcode))
            print('residues: {}'.format(self.atoms))
            print('frames: {}'.format(self.frames))
            print('temperature: {}'.format(self.temperature))

            pt.principal_axes(self.traj, dorotation=True)
            pt.center(self.traj, center='origin')

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

            # calculate fluctuations
            self.calcFluct()
            print('calculated fluctuations')
            # calculate correlations
            self.calcCorr()
            print('calculated correlations')
            # calculate principal axes (length)
            self.calcPrincipalAxes()
            print('calculated principal axes')
            # calculate radius of gyrations
            self.rg = np.mean(pt.radgyr(self.traj))
            # computes average distance between residues (~3.8A)
            self.avg_dist = np.mean([self.dist[i][i+1] for i in range(self.traj.n_atoms - 1)])
            # calculate shape factor fluctuations and mean value
            self.shape_factor_fluct()
            self.shape_factor = np.mean(self.SF_fluct)
            # compute suceptibility
            self.suceptibility = self.suceptibility()
            print('calculated suceptibility')

            if save_new:
                self.saveData()

    def saveData(self):
        f = open(self.data_path, 'wb')

        dataset = {'PDBcode': self.PDBcode,
                   'resnum': self.atoms,
                   'temperature': self.temperature,
                   'frames': self.frames,
                   'radgyr': self.rg,
                   'avg_dist': self.avg_dist,
                   'shape_factor': self.shape_factor,
                   'SF_fluct': self.SF_fluct,
                   'suceptibility': self.suceptibility,
                   'correl_len': self.corr_len,
                   'fluct_norms': self.fluct_norms,
                   'avgd_correlation': self.avgd_corr.to_dict(),
                   'correlation': self.corr,
                   'distance': self.dist,
                   'covariance': self.cov,
                   'fluctuations': self.fluct,
                   'trajectory': self.traj.xyz,
                   'ppal': self.ppal}

        pkl.dump(dataset, f)
        print('Saved data_{}'.format(self.PDBcode))
        f.close()

    def loadtraj(self):
        filesPath = os.path.join(self.prot, 'temp{}k'.format(self.temperature), 'input')
        trajPath = os.path.join(filesPath, 'nowat.{}.nc'.format(self.prot))
        topPath = os.path.join(filesPath, 'nowat.{}.prmtop'.format(self.prot))
        self.traj = pt.load(trajPath, topPath, mask='@CA', stride=self.strd)

    def avg_corr(self):
        flat_corr = np.ndarray.flatten(np.array(self.corr))
        flat_dist = np.ndarray.flatten(np.array(self.dist))

        self.avgd_corr = pd.DataFrame({'distance': flat_dist, 'correlation': flat_corr})
        self.avgd_corr.sort_values(by=['distance'], inplace=True)
        self.avgd_corr = self.avgd_corr[self.avgd_corr.distance != 0]
        self.avgd_corr = self.avgd_corr.rolling(self.atoms).mean()
        self.avgd_corr.dropna(axis=0, inplace=True)
        # correlation length defined as the first negative or zero value
        cl_idx = (self.avgd_corr.correlation <= 0).idxmax()
        self.corr_len = self.avgd_corr.distance[cl_idx]
        print('averaged correlation')

    def calcFluct(self):
        # First rmsd fit
        pt.rmsd(self.traj)
        # Calculate average structure
        self.avg = pt.mean_structure(self.traj)
        # Perform rms fit to average structure
        pt.rmsd(self.traj, ref=self.avg)
        # Calculate fluctuations
        for frame in self.traj:
            diff = frame.xyz - self.avg
            self.fluct.append(diff)
            self.fluct_norms.append([ln.norm(r) for r in diff])

    def calcCorr(self):
        # Compute covariance coefficient matrix
        for i in range(self.atoms):
            for j in range(self.atoms):
                for frm in self.fluct:
                    self.cov[i][j] += np.sum([frm[i][q] * frm[j][q] for q in range(3)])
                self.cov[i][j] /= self.frames

        # Compute correlation coefficient matrix
        # Calculate mutual distances between residues
        for i in range(self.atoms):
            for j in range(self.atoms):
                self.corr[i][j] = self.cov[i][j] / np.sqrt(self.cov[i][i] * self.cov[j][j])
                self.dist[i][j] = pt.distance(self.traj, ':{} :{}'.format(i+1, j+1), frame_indices=[0])[0]

    def calcPrincipalAxes(self):
        # Principal axes length
        xyz = self.traj.xyz

        minx = np.mean([min([q[0] for q in xyz[i]]) for i in range(self.frames)])
        maxx = np.mean([max([q[0] for q in xyz[i]]) for i in range(self.frames)])

        miny = np.mean([min([q[1] for q in xyz[i]]) for i in range(self.frames)])
        maxy = np.mean([max([q[1] for q in xyz[i]]) for i in range(self.frames)])

        minz = np.mean([min([q[2] for q in xyz[i]]) for i in range(self.frames)])
        maxz = np.mean([max([q[2] for q in xyz[i]]) for i in range(self.frames)])

        axisX = abs(maxx - minx) / 2
        axisY = abs(maxy - miny) / 2
        axisZ = abs(maxz - minz) / 2

        self.ppal = [axisX, axisY, axisZ]

    def shape_factor_fluct(self):
        xyz = self.traj.xyz

        for frm in xyz:
            Lx = abs(max([q[0] for q in frm]) - min([q[0] for q in frm])) / 2
            Ly = abs(max([q[1] for q in frm]) - min([q[1] for q in frm])) / 2
            Lz = abs(max([q[2] for q in frm]) - min([q[2] for q in frm])) / 2
            S = self.atoms * (3.8 ** 3) / (Lx * Ly * Lz)
            self.SF_fluct.append(S)

    def suceptibility(self):
        self.avg_corr()
        x = self.avgd_corr['distance'][self.avgd_corr['distance'] <= self.corr_len]
        y = self.avgd_corr['correlation'][self.avgd_corr['distance'] <= self.corr_len]
        suceptibility = (self.shape_factor / self.traj.n_atoms) * sum(y)
        return suceptibility
