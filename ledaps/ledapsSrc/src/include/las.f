#ifdef gould
	implicit undefined (a-z)
	include '$LASINCLUDE/worgen.fin'
	include '$LASINCLUDE/display.fin'    
	include '$LASINCLUDE/ddr.fin'
	include '$LASINCLUDE/imageiog.fin'
	include '$TAEINC/pgminc.fin'
#endif

#ifdef sun
	implicit undefined (a-z)
	include '$LASINCLUDE/worgen.fin'
	include '$LASINCLUDE/display.fin'    
	include '$LASINCLUDE/ddr.fin'
	include '$LASINCLUDE/imageios.fin'
	include '$TAEINC/pgminc.fin'
#endif

#ifdef vms
	implicit none
	include 'lasinclude:worgen.fin'
	include 'lasinclude:display.fin'
	include 'lasinclude:ddr.fin'
	include 'lasinclude:imageiov.fin'
	include 'tae$inc:pgminc.fin'
#endif
