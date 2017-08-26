```
# This file is part of the master thesis "Variational crimes in the Localized orthogonal decomposition method":
#   https://github.com/TiKeil/Masterthesis-LOD.git
# Copyright holder: Tim Keil 
# License: BSD 2-Clause License (http://opensource.org/licenses/BSD-2-Clause)
```

import numpy as np
import csv

def safeChange(ROOT, R, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest):
    C1Base = R.flatten()
    with open("%s/C1Base.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in C1Base:
            writer.writerow([val])

    with open("%s/C1vis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])

    #safe eps
    with open("%s/C1eps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])

    #safe PotentialUpdated
    with open("%s/C1PotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])

    #safe
    with open("%s/C1recomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])

    #safe
    with open("%s/C1errorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])

    #safe
    with open("%s/C1errorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])

    #safe
    with open("%s/C1errorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])
            
def safeVanish(ROOT, R, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest):
    VBase = R.flatten()
    with open("%s/VBase.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in VBase:
            writer.writerow([val])

    #safe vis
    with open("%s/Vvis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])

    #safe eps
    with open("%s/Veps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])

    #safe PotentialUpdated
    with open("%s/VPotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])

    #safe 
    with open("%s/Vrecomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])

    #safe 
    with open("%s/Verrorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])

    #safe 
    with open("%s/Verrorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])

    #safe 
    with open("%s/Verrorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])

def safeShift(ROOT, R, vis, eps, PotentialUpdated, recomputefractionsafe, errorplotinfo, errorworst, errorbest):
    M1Base = R.flatten()
    with open("%s/M1Base.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in M1Base:
            writer.writerow([val])


    with open("%s/M1vis.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in vis:
            writer.writerow([val])

    #safe eps
    with open("%s/M1eps.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in eps:
            writer.writerow([val])

    #safe PotentialUpdated
    with open("%s/M1PotentialUpdated.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in [PotentialUpdated]:
            writer.writerow([val])

    #safe 
    with open("%s/M1recomputefractionsafe.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in recomputefractionsafe:
            writer.writerow([val])

    #safe 
    with open("%s/M1errorplotinfo.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorplotinfo:
            writer.writerow([val])

    #safe 
    with open("%s/M1errorworst.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorworst:
            writer.writerow([val])

    #safe 
    with open("%s/M1errorbest.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in errorbest:
            writer.writerow([val])
            
def RegainChange(ROOT):
    C1Base = []
    f = open('%s/C1Base.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1Base.append(float(row[0]))
    f.close()

    C1eps = []
    f = open('%s/C1eps.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1eps.append(float(row[0]))
    f.close()

    C1errorbest = []
    f = open('%s/C1errorbest.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1errorbest.append(float(row[0]))
    f.close()

    C1errorplotinfo = []
    f = open('%s/C1errorplotinfo.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1errorplotinfo.append(float(row[0]))
    f.close()

    C1errorworst = []
    f = open('%s/C1errorworst.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1errorworst.append(float(row[0]))
    f.close()

    C1vis = []
    f = open('%s/C1vis.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1vis.append(float(row[0]))
    f.close()

    C1recomputefractionsafe = []
    f = open('%s/C1recomputefractionsafe.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        C1recomputefractionsafe.append(float(row[0]))
    f.close()
    
    return C1Base, C1eps, C1errorbest, C1errorplotinfo, C1errorworst, C1vis, C1recomputefractionsafe
    
def RegainVanish(ROOT):
    VBase = []
    f = open('%s/VBase.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        VBase.append(float(row[0]))
    f.close()

    Veps = []
    f = open('%s/Veps.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        Veps.append(float(row[0]))
    f.close()

    Verrorbest = []
    f = open('%s/Verrorbest.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        Verrorbest.append(float(row[0]))
    f.close()

    Verrorplotinfo = []
    f = open('%s/Verrorplotinfo.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        Verrorplotinfo.append(float(row[0]))
    f.close()

    Verrorworst = []
    f = open('%s/Verrorworst.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        Verrorworst.append(float(row[0]))
    f.close()

    Vvis = []
    f = open('%s/Vvis.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for rovw in reader:
        Vvis.append(float(row[0]))
    f.close()

    Vrecomputefractionsafe = []
    f = open('%s/Vrecomputefractionsafe.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        Vrecomputefractionsafe.append(float(row[0]))
    f.close()
    
    return VBase, Veps, Verrorbest, Verrorplotinfo, Verrorworst, Vvis, Vrecomputefractionsafe
    
def RegainShift(ROOT):
    M1Base = []
    f = open('%s/M1Base.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1Base.append(float(row[0]))
    f.close()

    M1eps = []
    f = open('%s/M1eps.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1eps.append(float(row[0]))
    f.close()

    M1errorbest = []
    f = open('%s/M1errorbest.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1errorbest.append(float(row[0]))
    f.close()

    M1errorplotinfo = []
    f = open('%s/M1errorplotinfo.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1errorplotinfo.append(float(row[0]))
    f.close()

    M1errorworst = []
    f = open('%s/M1errorworst.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1errorworst.append(float(row[0]))
    f.close()

    M1vis = []
    f = open('%s/M1vis.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1vis.append(float(row[0]))
    f.close()

    M1recomputefractionsafe = []
    f = open('%s/M1recomputefractionsafe.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        M1recomputefractionsafe.append(float(row[0]))
    f.close()
    
    return M1Base, M1eps, M1errorbest, M1errorplotinfo, M1errorworst, M1vis, M1recomputefractionsafe
    
def safer(ROOT, mum, a, plottingx, plottingy, plottingz, plotting2x, plotting2y, plotting2z, plotting3x, plotting3y, plotting3z, ems, Matrix):
    with open("%s/mum.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in mum:
            writer.writerow([val])
    
    with open("%s/a.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in a:
            writer.writerow([val])
    
    with open("%s/Matrix.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in Matrix:
            writer.writerow([val])
            
    with open("%s/ems.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in ems:
            writer.writerow([val])
    
    with open("%s/plottingx.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plottingx:
            writer.writerow([val])

    with open("%s/plottingy.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plottingy:
            writer.writerow([val])
    
    with open("%s/plottingz.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plottingz:
            writer.writerow([val])
            
    with open("%s/plotting2x.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting2x:
            writer.writerow([val])

    with open("%s/plotting2y.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting2y:
            writer.writerow([val])
    
    with open("%s/plotting2z.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting2z:
            writer.writerow([val])
    
    with open("%s/plotting3x.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting3x:
            writer.writerow([val])

    with open("%s/plotting3y.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting3y:
            writer.writerow([val])
    
    with open("%s/plotting3z.txt" % ROOT, 'wb') as csvfile:
        writer = csv.writer(csvfile)
        for val in plotting3z:
            writer.writerow([val])
    
def regainer(ROOT):
    mum = []
    f = open('%s/mum.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        mum.append(float(row[0]))
    f.close()
    
    a = []
    f = open('%s/a.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        a.append(float(row[0]))
    f.close()

    ABase = []
    f = open('%s/Matrix.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        ABase.append(float(row[0]))
    f.close()
    
    ems = []
    f = open('%s/ems.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        ems.append(float(row[0]))
    f.close()
    
    plottingx = []
    f = open('%s/plottingx.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plottingx.append(float(row[0]))
    f.close()

    plottingy = []
    f = open('%s/plottingy.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plottingy.append(float(row[0]))
    f.close()
    
    plottingz = []
    f = open('%s/plottingz.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plottingz.append(float(row[0]))
    f.close()

    plotting2x = []
    f = open('%s/plotting2x.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting2x.append(float(row[0]))
    f.close()

    plotting2y = []
    f = open('%s/plotting2y.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting2y.append(float(row[0]))
    f.close()
    
    plotting2z = []
    f = open('%s/plotting2z.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting2z.append(float(row[0]))
    f.close()
    
    plotting3x = []
    f = open('%s/plotting3x.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting3x.append(float(row[0]))
    f.close()

    plotting3y = []
    f = open('%s/plotting3y.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting3y.append(float(row[0]))
    f.close()
    
    plotting3z = []
    f = open('%s/plotting3z.txt' % ROOT, 'rb')
    reader = csv.reader(f)
    for row in reader:
        plotting3z.append(float(row[0]))
    f.close()
    
    return mum, a, ABase, ems, plottingx, plottingy, plottingz, plotting2x, plotting2y, plotting2z, plotting3x, plotting3y, plotting3z    