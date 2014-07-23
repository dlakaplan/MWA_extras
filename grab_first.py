datacol='data'

ms.open('1063400456_corr.ms')
scansummary = ms.getscansummary()
starttime_mjd = scansummary['1']['0']['BeginTime']
starttime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0/(24.*60*60),'d'),form=['ymd'], prec=9)[0], 's'))[0]
stoptime0 = qa.getvalue(qa.convert(qa.time(qa.quantity(starttime_mjd+0.5/(24.*60*60), 'd'), form=['ymd'], prec=9)[0], 's'))[0]

ms.selectinit(datadescid=0)  # initialize to initialize params
selection = {'time': [starttime0, stoptime0]}
ms.select(items = selection)
da = ms.getdata([datacol, 'axis_info','u','v','w'], ifraxis=True)
ms.close()
