/*globals Peak, nmr*/
function printPeakAttrs(){
    
    print("Peak attributes:");
    print("Peak.Type.Artifact " + Peak.Types.Artifact );
    print("Peak.Types.CS_Ref " + Peak.Types.CS_Ref );
    print("Peak.Types.Compound " + Peak.Types.Compound );
    print("Peak.Types.Impurity " + Peak.Types.Impurity );
    print("Peak.Types.Quantification_Ref " + Peak.Types.Quantification_Ref );
    print("Peak.Types.Solvent " + Peak.Types.Solvent );

    print("\nPeak Flags\n" );

    // C13Satellite, Hidden, Labile, None, Rotational, Solvent2ary, Weak
    print("Peak.Flags.C13Satellite " + Peak.Flags.C13Satellite );
    print("Peak.Flags.Hidden " + Peak.Flags.Hidden );
    print("Peak.Flags.Labile " + Peak.Flags.Labile );
    print("Peak.Flags.None " + Peak.Flags.None );
    print("Peak.Flags.Rotational " + Peak.Flags.Rotational );
    print("Peak.Flags.Solvent2ary " + Peak.Flags.Solvent2ary );
    print("Peak.Flags.Weak " + Peak.Flags.Weak );

    print("\nPeak Kinds\n" );

    // Maximum, Minimum
    print("Peak.Kinds.Maximum " + Peak.Kinds.Maximum );
    print("Peak.Kinds.Minimum " + Peak.Kinds.Minimum );

}