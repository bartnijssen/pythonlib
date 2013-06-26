import datetime as dt
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import vic.IO
from utils.plot import cmap_discretize


def parsepolyfile(polyfile):
    '''Read file with polynomial coefficients'''
    with open(polyfile, 'r') as f:
        order = int(f.readline().split()[0])
        polyparams = []
        for i in range(order+1):
            polyparams.append(float(f.readline()))
        fields = f.readline().split()
        (xmin, xmax) = [float(i) for i in fields[:2]]
    return (np.poly1d(polyparams), xmin, xmax)


def determine_update_function(station, timeofconcentration, startwindow=dt.datetime(1950, 1, 1),
                              endwindow=dt.datetime(2100, 12, 31),
                              masktemplate='/Users/nijssen/data/pnnl/test_new_direction/get_maskfile/{}.ll',
                              fluxtemplate='/Users/nijssen/data/pnnl/19492010/fluxes/fluxes_{}_{}',
                              vicflowtemplate='/Users/nijssen/data/pnnl/19492010/flow_2/{}.day',
                              polytemplate='/Users/nijssen/data/pnnl/poly/{}.poly',
                              figuretemplate=None,
                              latlondigits=4):
    '''Determine the polynomial fit between soil moisture and runoff'''
    rollingwindow = timeofconcentration + 1
    # Read the flux files for all cells in the basin
    fluxtemplate = fluxtemplate.format('{{:.{}f}}'.format(latlondigits),
                                       '{{:.{}f}}'.format(latlondigits))
    maskfile = masktemplate.format(station)
    with open(maskfile, 'r') as f:
        read_data = f.readlines()
    ncells = 0
    for line in read_data:
        (lon, lat) = [float(i) for i in line.split()]
        nf = vic.IO.readfluxfile(fluxtemplate.format(lat, lon))
        if ncells == 0:
            vicsm = nf
        else:
            vicsm += nf
        ncells += 1
    vicsm /= ncells
    vicsm = pd.rolling_mean(vicsm, rollingwindow)
    vicsm = vicsm.tshift(-rollingwindow, freq='D')
    vicsm = vicsm[startwindow:]

    # Read the simulated (and routed) flow data
    flowfile = vicflowtemplate.format(station)
    vicflow = vic.IO.readflowfile(flowfile)
    vicflow = pd.rolling_mean(vicflow, rollingwindow)
    vicflow = vicflow.tshift(-rollingwindow, freq='D')
    vicflow = vicflow[startwindow:]

    # Determine the third order polynomial relating SM_3 and R
    idx = ((~np.isnan(vicsm['sm3'])) & (~np.isinf(vicsm['sm3'])) &
           (~np.isnan(vicflow['flow'])) & (~np.isinf(vicflow['flow'])))
    Y = vicflow['flow'][idx]
    X = vicsm['sm3'][idx]
    poly = np.poly1d(np.polyfit(X, Y, 3))

    # Write the polynomial to file
    polyfile = polytemplate.format(station)
    with open(polyfile, 'w') as f:
        f.write('{:d} # order of polynomial for station: '.format(poly.order))
        f.write('{} with {}-day time of concentration\n'.format(station, timeofconcentration))
        np.savetxt(f, poly)
        f.write('{} {} # range of soil moisture values\n'.format(X.min(), X.max()))
        f.write('###########################\n')
        f.write('Polynomial in readable form:\n')
        f.write(str(poly))

    if figuretemplate is not None:
        # Create a diagnostic plot of the relationship between SM_3 and R
        figurefile = figuretemplate.format(station)
        colormap = cmap_discretize(mpl.cm.get_cmap('RdYlBu'), 12)
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(111)
        plt.scatter(X, Y, marker='.', c=X.index.month, edgecolor='none',
                    cmap=colormap)
        cbar = plt.colorbar()
        cbar.set_label('Month')
        cbar.set_ticks((np.arange(0, 13)+0.5)*12/13)
        cbar.set_ticklabels(['{:02d}'.format(i) for i in range(1, 13)])
        plt.scatter(X, poly(X), marker='.', edgecolor='none', c='b')
        plt.ylim(ymin=-100)
        plt.xlabel('Soil moisture layer 3 (mm)')
        plt.ylabel('Simulated flow (cfs)')
        plt.grid(True)
        ax.text(0.03, 0.97, station + ', window size: {} days'.format(rollingwindow),
                verticalalignment='top', horizontalalignment='left',
                transform=ax.transAxes, color='black', fontsize=16)
        polystr = '$Q='
        for i in range(poly.order):
            polystr += '{:.4g}'.format(poly.c[i])
            if (poly.order-i) > 1:
                polystr += 'SM_3^{}+'.format(poly.order-i)
            else:
                polystr += 'SM_3+'
        polystr += '{:.4g}$'.format(poly.c[poly.order])
        ax.text(0.03, 0.90, r'{}'.format(polystr), verticalalignment='top', horizontalalignment='left',
                transform=ax.transAxes, color='black', fontsize=14)

        plt.savefig(figurefile)

    return poly


def assimilate_streamflow(forecastdate=dt.datetime(2012, 4, 21),
                          stationdict={'NFORK': 2, 'UNCOM': 2},
                          statetemplate='/Users/nijssen/data/pnnl/19492010/states/state_{:04d}{:02d}{:02d}',
                          soilfile='/Users/nijssen/data/pnnl/projects/waterpower/forecast/uswide/data/colo/params/vic/calibrated_10_04.soil.fix.fc.colo',
                          masktemplate='/Users/nijssen/data/pnnl/test_new_direction/get_maskfile/{}.ll',
                          fluxtemplate='/Users/nijssen/data/pnnl/19492010/20092010fluxes/fluxes_{:.4f}_{:.4f}',
                          vicflowtemplate='/Users/nijssen/data/pnnl/19492010/flow_2/{}.day',
                          obsflowtemplate='/Users/nijssen/data/pnnl/19492010/flow_obs/{}.day.historic',
                          polytemplate='/Users/nijssen/data/pnnl/poly/{}.poly',
                          latlondigits=4):
    '''Update soil moisture based on observed streamflow'''
    # As explained above, I will only handle cases for which the time of
    # concentration is the same for all basins (since that means I only
    # need to touch a single state file

    timeofconcentration = None
    for station, tc in stationdict.iteritems():
        if timeofconcentration is None:
            timeofconcentration = tc
        else:
            if tc != timeofconcentration:
                raise ValueError("Time of concentration expected to be the same for all stations")
    rollingwindow = timeofconcentration + 1

    # Determine the date of the update. This will be tc + 1 days before the forecast
    updatedate = forecastdate - dt.timedelta(days=rollingwindow)

    # For each station create a list with grid cells associated with that station and
    # for each grid cell record the associated station (to facilitate two-way lookup)
    cellstationdict = {}
    stationinfo = {}
    celltemplate = '{{:.{}f}}_{{:.{}f}}'.format(latlondigits, latlondigits)
    for station in stationdict.iterkeys():
        stationinfo[station] = {}
        stationinfo[station]['cells'] = []
        maskfile = masktemplate.format(station)
        with open(maskfile, 'r') as f:
            read_data = f.readlines()
        for line in read_data:
            (lon, lat) = [float(i) for i in line.split()]
            cell = celltemplate.format(lat, lon)
            if cell not in cellstationdict:
                cellstationdict[cell] = []
            cellstationdict[cell].append(station)
            stationinfo[station]['cells'].append(cell)

    # For each station parse the poly file and store as poly1d object
    for station in stationdict.iterkeys():
        (stationinfo[station]['poly'],
         stationinfo[station]['minsm'],
         stationinfo[station]['maxsm']) = parsepolyfile(polytemplate.format(station))
        stationinfo[station]['minflow'] = stationinfo[station]['poly'](stationinfo[station]['minsm'])
        stationinfo[station]['maxflow'] = stationinfo[station]['poly'](stationinfo[station]['maxsm'])

    # For each station parse the flux files and store the mean soil moisture in the third layer
    fluxtemplate = fluxtemplate.format('{{:.{}f}}'.format(latlondigits),
                                       '{{:.{}f}}'.format(latlondigits))
    for station in stationdict.iterkeys():
        ncells = 0
        for cell in stationinfo[station]['cells']:
            (lat, lon) = [float(i) for i in cell.split('_')]
            nf = vic.IO.readfluxfile(fluxtemplate.format(lat, lon))
            if ncells == 0:
                vicsm = nf
            else:
                vicsm += nf
            ncells += 1
        vicsm /= ncells
        vicsm = pd.rolling_mean(vicsm, rollingwindow)
        # time shift the window, so it is forward-looking instead of backward-looking
        vicsm = vicsm.tshift(-rollingwindow, freq='D')
        stationinfo[station]['oldsm'] = vicsm['sm3'][updatedate]

    # Parse the soil file and store the minimum and maximum value for SM3
    cellsm3limits = vic.IO.getsmlimits(soilfile, cellstationdict.keys(), celltemplate)

    # Parse the state file for the update date
    state = vic.IO.parsestatefile(statetemplate.format(updatedate.year, updatedate.month, updatedate.day))

    # Parse the routed and observed flows
    for station in stationdict.iterkeys():
        vicflow = vic.IO.readflowfile(vicflowtemplate.format(station))
        vicflow = pd.rolling_mean(vicflow, rollingwindow)
        vicflow = vicflow.tshift(-rollingwindow, freq='D')
        stationinfo[station]['vicflow'] = vicflow['flow'][updatedate]
        obsflow = vic.IO.readflowfile(obsflowtemplate.format(station), skiprows=1, conversion=(1/.3048)**3)
        obsflow = pd.rolling_mean(obsflow, rollingwindow)
        obsflow = obsflow.tshift(-rollingwindow, freq='D')
        stationinfo[station]['obsflow'] = obsflow['flow'][updatedate]

    # all the relevant info has been read. Now determine the update
    for station in stationdict.iterkeys():
        # If the observed flow is greater than the maximum flow with the polynomial, set the new soil moisture
        # to the maximum of the modeled soil moisture and the range used in determining the polynomial
        if stationinfo[station]['obsflow'] > stationinfo[station]['maxflow']:
            stationinfo[station]['newsm'] = max(stationinfo[station]['oldsm'], stationinfo[station]['maxsm'])
        # If the observed flow is less than the minimum flow with the polynomial, set the new soil moisture
        # to the minimum of the modeled soil moisture and the range used in determining the polynomial
        elif stationinfo[station]['obsflow'] < stationinfo[station]['minflow']:
            stationinfo[station]['newsm'] = min(stationinfo[station]['oldsm'], stationinfo[station]['minsm'])
        # find the roots of the polynomial and select the one that is within the range
        else:
            updatepoly = stationinfo[station]['poly'] - stationinfo[station]['obsflow']
            roots = updatepoly.roots
            newsm = np.real(roots[(roots >= stationinfo[station]['minsm']) &
                                  (roots <= stationinfo[station]['maxsm']) &
                                  (np.isreal(roots))])[0]
            # do NOT update if the update takes you in the wrong direction
            if ((stationinfo[station]['obsflow'] < stationinfo[station]['vicflow']) &
                (newsm > stationinfo[station]['oldsm'])):
                # do nothing, which means set the new soil to the old one
                stationinfo[station]['newsm'] = stationinfo[station]['oldsm']
            elif ((stationinfo[station]['obsflow'] > stationinfo[station]['vicflow']) &
                  (newsm < stationinfo[station]['oldsm'])):
                stationinfo[station]['newsm'] = stationinfo[station]['oldsm']
            else:
                stationinfo[station]['newsm'] = newsm

        stationinfo[station]['updateratio'] = stationinfo[station]['newsm'] / stationinfo[station]['oldsm']

    # update the state
    for cell in cellstationdict.iterkeys():
        cellid = cellsm3limits[cell][0]
        updateratio = 0
        # if there are multiple stations associated with the same cell, average the update
        for station in cellstationdict[cell]:
            updateratio += stationinfo[station]['updateratio']
        updateratio /= len(cellstationdict[cell])
        if updateratio != 1:
            nveg = state['cells'][cellid]['nveg']
            nbands = state['cells'][cellid]['nbands']
            for i in range(nveg+1):
                for j in range(nbands):
                    oldsm = float(state['cells'][cellid]['vegband_data'][i][j][4])
                    newsm = updateratio * oldsm
                    # it looks like the soil moisture in the version of VIC that is being used
                    # can fall below the minimum value, so even though I think the following
                    # line is correct, I'll replace the test with one for newsm > 0 for now
                    # newsm = max(newsm, cellsm3limits[cell][1])
                    newsm = max(newsm, 0)
                    newsm = min(newsm, cellsm3limits[cell][2])
                    state['cells'][cellid]['vegband_data'][i][j][4] = newsm

    # write the new state to file
    vic.IO.writestatefile(statetemplate.format(updatedate.year, updatedate.month, updatedate.day)+'.update', state)

    return stationinfo
