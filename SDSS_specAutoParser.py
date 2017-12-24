# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 17:40:44 2017

@author: Qing Liu
@institute: USTC

This program is a demo to auto-scrape the H_alpha line fluxes of 20 
randomly selected emission-line galaxies from the Line Measurement 
Information table listed on the SDSS Optical Spectra webpage, given 
PlateID, MJD and FiberID of each target spectrum.

(To include more lines, modify the beginning & end of the func)

"""

import numpy as np
import pandas as pd
import requests
from astropy.io import ascii
from astropy.table import Table
from bs4 import BeautifulSoup



def AutoScrape_SDSS_spline(Plateid, MJD, Fiberid, N):

	ha_flx = np.zeros(N) 	# set Ha

	for i, (mjd, fiber, plateid) in enumerate(zip(MJD, Fiberid, Plateid)):

		print("Requesting SDSS spectrum #%d... mjd=%s fiber=%s plateid=%s" %(i+1, mjd,fiber,plateid))

		url = 'https://dr10.sdss.org/spectrumDetail?mjd=%s&fiber=%s&plateid=%s' %(mjd,fiber,plateid)

		# request the url and return in html 
		r = requests.get(url)
		content = r.text

		# parse with BeautifulSoup & search pattern
		soup = BeautifulSoup(content, 'lxml')
		tabular = soup.find("table", 
				{"class": "table table-striped table-bordered", "id":"lines"})

		# (data start from 2nd row)
		tab_lines = tabular.find_all('tr')[1:] 

		# extract emission lines info
		line_flx = np.zeros(len(tab_lines))

		for j, tab_line in enumerate(tab_lines):
			line_info = tab_line.find_all("td")

			# column 1: line name
			line_name = line_info[0].get_text()

			# column 5: line flux
			line_flx[j] = float(line_info[4].get_text())
	
			if line_name == "H_alpha":
				ha_flx[i] = line_flx[j]

	return ha_flx


#####-----Program Start Here-----#####

# read global table
table_lib = ascii.read("./astro-lab.lis").to_pandas()

# random sampling
N = 20
table = table_lib.sample(N)

# run auto-scraping
Ha_flux = AutoScrape_SDSS_spline(Plateid = table.plateid, 
				 MJD = table.mjd, 
				 Fiberid = table.fiberid, N=N)

print("Finish Auto-scraping! Save to txt File")

# write data to txt file
data = Table([table.plateid, table.mjd, table.fiberid, Ha_flux],
		 names=['plateid', 'MJD', 'fiberid', 'Ha_flux'])

ascii.write(data, 'SDSS_LineInfo.txt', overwrite=True)	
