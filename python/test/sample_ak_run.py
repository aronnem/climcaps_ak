import numpy as np
import h5py

from pyclimcaps import return_climcaps_akdata

def write_ak_run_output(
        h5_output_file,
        Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel):
    with h5py.File(h5_output_file, 'w') as h:
        h['Fmatrix'] = Fmatrix
        h['Finv'] = Finv
        h['AKcoarse'] = AKcoarse
        h['Pcoarse'] = Pcoarse
        h['AKfine'] = AKfine
        h['Pfine'] = Pfine
        h['Skernel'] = Skernel


def test_ak_run_output(
        h5_output_file,
        Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel):
    with h5py.File(h5_output_file, 'r') as h:
        print('Fmatrix max difference:  {:15.7g}'.format(np.max(np.abs(Fmatrix - h['Fmatrix'][:]))))
        print('Finv max difference:     {:15.7g}'.format(np.max(np.abs(Finv - h['Finv'][:]))))
        print('AKcoarse max difference: {:15.7g}'.format(np.max(np.abs(AKcoarse - h['AKcoarse'][:]))))
        print('Pcoarse max difference:  {:15.7g}'.format(np.max(np.abs(Pcoarse - h['Pcoarse'][:]))))
        print('AKfine max difference:   {:15.7g}'.format(np.max(np.abs(AKfine - h['AKfine'][:]))))
        print('Pfine max difference:    {:15.7g}'.format(np.max(np.abs(Pfine - h['Pfine'][:]))))
        print('Skernel max difference:  {:15.7g}'.format(np.max(np.abs(Skernel - h['Skernel'][:]))))


def sample_ak_run(molnames, write_h5_output=False, test_h5_output=True):

    # for each of these two files, there is a hand-picked FOR.
    # these two files should be downloaded with the data_download.sh file
    # at the top level.

    ncfiles = [
        '../../test_data/SNDR.J1.CRIMSS.20190901T0336.m06.g037.L2_CLIMCAPS_RET.std.v02_28.G.200214174949.nc',
        '../../test_data/SNDR.J1.CRIMSS.20190901T2248.m06.g229.L2_CLIMCAPS_RET.std.v02_28.G.200214190447.nc' ]
    ifoot = [3, 4]
    iscan = [14, 32]

    for m,i in np.ndindex(len(molnames), len(ncfiles)):

        Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel = \
            return_climcaps_akdata(ncfiles[i], ifoot[i], iscan[i], molnames[m])

        #outname = f'for{ifoot[i]:d}_scan{iscan[i]:d}_{molnames[m]:s}'

        #plot_fmatrix, Fmatrix, Finv, Pfine, 'example_plot_Fmatrix_'+outname+'.png'
        #plot_ak, AKcoarse, Pcoarse, 'example_plot_AKcoarse_'+outname+'.png'
        #plot_ak, AKfine, Pfine, 'example_plot_AKfine_'+outname+'.png'
        #plot_sfunc, Skernel, Pfine, 'example_plot_Sfunc_'+outname+'.png'

        #write_csv, 'example_Fmatrix_'+outname+'.csv', Fmatrix
        #write_csv, 'example_Finv_'+outname+'.csv', Finv
        #write_csv, 'example_AKfine_'+outname+'.csv', AKfine
        #write_csv, 'example_AKcoarse_'+outname+'.csv', AKcoarse
        #write_csv, 'example_Skernel_'+outname+'.csv', Skernel

        h5_output_file = f'../../test_data/py_test_case_{i+1:1d}_{molnames[m]:s}.h5'
        h5_ref_file = f'../../test_data/test_case_{i+1:1d}_{molnames[m]:s}.h5'
        if write_h5_output:
            print('writing output file: ', h5_output_file)
            write_ak_run_output(
                h5_output_file, 
                Fmatrix, Finv, AKcoarse, Pcoarse, AKfine, Pfine, Skernel)
        if test_h5_output:
            print('testing output relative to: ', h5_ref_file)
            # reference data came from IDL, and has transposed arrays.
            test_ak_run_output(
                h5_ref_file,
                Fmatrix.T, Finv.T, AKcoarse.T, Pcoarse.T, AKfine.T, Pfine.T, Skernel.T)
