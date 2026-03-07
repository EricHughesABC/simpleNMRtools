function read_1D_eeh() {

  var spec = nmr.activeSpectrum();
  
  	if (!spec.isValid()) {
		return;
	}
	
	
	// print("dimCount " + spec.dimCount );	
	// print("spec.count " + spec.count(1) + " " + spec.count(2) );
	
	var rrr = spec.real({from:780, to:780}, {from:0,to:512});
	
	print( "rrr = [" );
	for( var i=0; i<spec.count(2); i++ ){
		print( rrr[i] + "," );
	}
	
	print( "]" );
	
	print(spec.real({from:0, to:10}, {from:0,to:10}));
	

	
	
	
			
}