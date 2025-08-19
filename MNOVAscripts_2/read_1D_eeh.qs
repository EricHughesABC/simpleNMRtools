function read_1D_eeh() {

  var spec = nmr.activeSpectrum();
  
  	if (!spec.isValid()) {
		return;
	}
	
	
	var sLimits = spec.getScaleLimits();
    	
	
	print("dimCount " + spec.dimCount );	
	print("spec.count " + spec.count(1) + " " + spec.count(2) );
	print("spec Threshold " + spec.threshold );
	print(" spec realMax " + spec.realMax );
	print("spec realMin " + spec.realMin );

	// print out the the sLimits
	//  sLimits.fromX, sLimits.toX, sLimits.fromY, sLimits.toY
	print("sLimits " + sLimits.fromX + " " + sLimits.toX + " " + sLimits.fromY + " " + sLimits.toY );
	
	
	print(sLimits);
	var rrr = spec.real({from:780, to:780}, {from:0,to:512});
	
	//print( "rrr = [" );
	//for( var i=0; i<512; i++ ){
	//	print( rrr[i] + "," );
	//}
	
	//print( "]" );
	
	//print(spec.real({from:0, to:10}, {from:0,to:10}));
	

	
	
	
			
}