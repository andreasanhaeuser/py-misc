#!/usr/bin/python3

import meteo_density as md

substance = 'O3'
mass_ug = 1

ppb = md.ug_to_ppb(mass_ug, substance)
print('ppb = ' + str(ppb))
