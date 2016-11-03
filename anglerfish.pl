#!/usr/bin/perl
use strict;
use List::Util qw( min max sum);
use Math::Trig;
use CGI;
my $head = qq{<head>
 <style>
  h1 {color: #46709D;}
  header {padding: 10Px;
    text-indent: 40px;
    max-width: 700px;
    font-family: "Lucida sans";}
  p {
    font-family: "Lucida sans";
    font-size: 10pt;
    max-width: 700px;
    }
  ul {font-family: "Lucida sans";
  font-size: 10pt;
    max-width: 750px;}
  table {font-family: "Lucida sans";}
  footer {color: gray;
    font: 10pt "Times New Roman", Times, serif;
    padding:20px;
    text-align:center;}
  </style>
 </head>
 <header>
  <h1><a href="../index.html"><img src="../logo.png" style="width:200px;height:128px;"></a>Anglerfish</h1>
   <p><center><a href="../help.html">Help</a> / <a href="../index.html">Main</a></center></p>
 </header>};
my $foot = qq{<footer> 
    <p> Anglerfish is designed by Matthew Colledge of the Wallace group<br>Institute of Structural and Molecular Biology, Birkbeck College, University of london</p>
  </footer>};
my $in = new CGI;
my $Ni = $in->param("Ni"); 
my $Ci = $in->param("Ci");
my $filename = $in->param("PDBfile");
my $warnings = "<b>Warnings and Errors</b><br>";
if ($Ci - $Ni <= 7) {
	$warnings .= "Input helix is less than 8 residues long, this calculation might not be sensical.<br>";
	} #checking numbers make sense 
if ($Ci - $Ni <= 0) {
	$warnings .= "C end of helix before N, helix must be specified in the N to C direction.<br>";
	}
my %x;
my %y;
my %z;
my $nextmodel = 2; #default of using model 1
my @chains;
my @residues;
my @lines = <$filename>;
foreach (@lines){ #find model no and chains
	if ($_ =~ /^REMARK 210 BEST REPRESENTATIVE CONFORMER IN THIS ENSEMBLE : (\d+)/){ #changes to best model
		$nextmodel = $1+1;
	}
	if ($_ =~ /^ATOM\s+\d+\s+CA\s+\w+\s+([A-Z])\s+(\d+)/) {
		push(@chains, $1) unless grep{$_ eq $1} @chains;
		push(@residues, $2) unless grep{$_ eq $2} @chains;
	}
}
if ($in->param("Modelno")) {
	$nextmodel = $in->param("Modelno")+1;
}
if ($in->param("Chains")) {
	if ($in->param("Chains") =~ /^([A-Z]\,){2,}[A-Z]\,?$/) {
		undef @chains;
		@chains = split (/,/,$in->param("Chains"));
	}
	else {
		$warnings .= "Chains to use entered incorretly, only type individual uppercase characters separated by single commas eg 'A,B,C' if you wish to use all the chains in the file you may leave this field empty.<br>";
		}
}
if ($in->param("Axi")) {
	@residues = ($in->param("Axi"));
}
foreach (@lines){
	if ($_ =~ /^ATOM\s+\d+\s+CA\s+\w+\s+(\w)\s*(\d+)\s+(-?\d+\.\d{3})\s*(-?\d+\.\d{3})\s*(-?\d+\.\d{3})/){ 
        	$x{$1,$2} = $3;
        	$y{$1,$2} = $4;
        	$z{$1,$2} = $5;
        	}
	if ($_ =~ /^MODEL\w+$nextmodel/){
    		last;
    		}
	}
foreach (@chains) {
	unless ($x{$_,$Ni}) {$warnings .= "Start of helix (residue $Ni) not found in chain $_<br>";}
	unless ($x{$_,$Ci}) {$warnings .= "End of helix (residue $Ci) not found in chain $_<br>";}
	}
unless (%x) {
	$warnings .= "No lines matching PDB format found in file, are you sure this is the correct file?<br>";
	print <<ENDERROR;
Content-type: text/html

	
<!DOCTYPE HTML>
<html>
$head
 <body>
 <p>$warnings</p>
 </body>
$foot
</html>
ENDERROR
	die;
}
my $model = $nextmodel-1;
my %RfA = ( x => 1, y => 1, z => 1); #reference cordinates for ABC to define up
$RfA{x} = $x{"$chains[0]",$residues[0]};
$RfA{y} = $y{"$chains[0]",$residues[0]};
$RfA{z} = $z{"$chains[0]",$residues[0]};
my %RfB; ( x => 1, y => 1, z => 1);
$RfB{x} = $x{"$chains[1]",$residues[0]};
$RfB{y} = $y{"$chains[1]",$residues[0]};
$RfB{z} = $z{"$chains[1]",$residues[0]};
my %RfC; ( x => 1, y => 1, z => 1);
$RfC{x} = $x{"$chains[2]",$residues[0]};
$RfC{y} = $y{"$chains[2]",$residues[0]};
$RfC{z} = $z{"$chains[2]",$residues[0]};
my %RfAtoB; # vectors from A to B and C defined to calculate the normal to the plane containing (the N terminal ends of) A, B and C which I hope is parralel to the is for a symmetric pore 
$RfAtoB{x} = $RfB{x} - $RfA{x};
$RfAtoB{y} = $RfB{y} - $RfA{y};
$RfAtoB{z} = $RfB{z} - $RfA{z};
my %RfAtoC;
$RfAtoC{x} = $RfC{x} - $RfA{x};
$RfAtoC{y} = $RfC{y} - $RfA{y};
$RfAtoC{z} = $RfC{z} - $RfA{z};
my %RfAtoAxis; #is defined as the normal vector at the Axi on chain A to the plane made by vectors to Axi at chains B and C
$RfAtoAxis{x} = $RfAtoB{y} * $RfAtoC{z} - $RfAtoB{z} * $RfAtoC{y};
$RfAtoAxis{y} = $RfAtoB{z} * $RfAtoC{x} - $RfAtoB{x} * $RfAtoC{z};
$RfAtoAxis{z} = $RfAtoB{x} * $RfAtoC{y} - $RfAtoB{y} * $RfAtoC{x};
my $RfAtoAxismag = sqrt(($RfAtoAxis{x}**2)+($RfAtoAxis{y}**2)+($RfAtoAxis{z}**2)); #reference vector
my @Axchains = @chains;
if ($in->param("Axchains")) {
	if ($in->param("Axchains") =~ /^([A-Z]\,){2,}[A-Z]\,?$/) {
		undef @Axchains;
		@Axchains = split (/,/,$in->param("Axchains"));
	}
	else {
		$warnings .= "Chains to define axis entered incorretly, only type individual uppercase characters separated by single commas eg 'A,B,C' if you wish to use all the chains in the file you may leave this field empty.<br>";
	}
} 
my $Axmer = scalar @Axchains; #number of chains to define axis
my $mer = scalar @chains; #mer as in order of homomer (number of chains)
unless ($Axmer > 2) {$warnings .= "Less than 3 chains found in PDB file, atleast 3 must be present and symmetric for the calculation.<br>";} 
my @Axisx; #arrays for the chosen axis average coordinates
my @Axisy;
my @Axisz;
my @Axangles; #array to check angle of axis to plane (should be 90 or -90)
my @notAxangles; #array to check discarded (should be 90 or -90)
my $dev = 1e99;
my $selaxi;
foreach my $Axi (@residues) {
	my @tAxisx=();
	my @tAxisy=();
	my @tAxisz=();
	for (my $i=0; $i < ($Axmer); $i++) {
		my %AxA; #store axis residue coordinates now just for check
		my %AxB;
		my %AxC;
		$AxA{x} = $x{$Axchains[$i],$Axi};
		$AxA{y} = $y{$Axchains[$i],$Axi};
		$AxA{z} = $z{$Axchains[$i],$Axi};
		unless ($AxA{x} and $AxA{y} and $AxA{z}){
							$warnings .= qq{Failed to find selected axis residue ($Axi) in chain $Axchains[$i], this residue must be present in all chains selected and distributed evenly around the axis of rotation for the calculation.<br>To omit chains from the calculation enter only the cahins you would like to use in the "chains to use" field.<br>};
						}
		for (my $j=0; $j < ($Axmer); $j++){
			unless ($i == $j) {
				$AxB{x} = $x{$Axchains[$j],$Axi};
				$AxB{y} = $y{$Axchains[$j],$Axi};
				$AxB{z} = $z{$Axchains[$j],$Axi};
				for (my $k=0; $k < ($Axmer); $k++){
					unless ($i == $k or $j == $k){
						$AxC{x} = $x{$Axchains[$k],$Axi};
						$AxC{y} = $y{$Axchains[$k],$Axi};
						$AxC{z} = $z{$Axchains[$k],$Axi};
						my %AtoB; # vectors from A to B and C defined to calculate the normal to the plane containing (the N terminal ends of) A, B and C which I hope is parralel to the is for a symmetric pore 
						$AtoB{x} = $AxB{x} - $AxA{x};
						$AtoB{y} = $AxB{y} - $AxA{y};
						$AtoB{z} = $AxB{z} - $AxA{z};
						my %AtoC;
						$AtoC{x} = $AxC{x} - $AxA{x};
						$AtoC{y} = $AxC{y} - $AxA{y};
						$AtoC{z} = $AxC{z} - $AxA{z};
						my %AxAtoAxis; #is defined as the normal vector at the Axi on chain A to the plane made by vectors to Axi at chains B and C
						$AxAtoAxis{x} = $AtoB{y} * $AtoC{z} - $AtoB{z} * $AtoC{y};
						$AxAtoAxis{y} = $AtoB{z} * $AtoC{x} - $AtoB{x} * $AtoC{z};
						$AxAtoAxis{z} = $AtoB{x} * $AtoC{y} - $AtoB{y} * $AtoC{x};
						unless ($AxAtoAxis{x} or $AxAtoAxis{y} or $AxAtoAxis{z}){
							print <<ENDERROR;
Content-type: text/html

	
<!DOCTYPE HTML>
<html>
$head
<body><p>$warnings</p>
</body>
$foot
</html>
ENDERROR
							die;
						}
						my $AxAtoAxismag = sqrt(($AxAtoAxis{x}**2)+($AxAtoAxis{y}**2)+($AxAtoAxis{z}**2)); #calculating the angle between the line AxAtoAxC and AxA to the new axis at AxA
						my $dotproduct = $AxAtoAxis{x}*$RfAtoAxis{x}+$AxAtoAxis{y}*$RfAtoAxis{y}+$AxAtoAxis{z}*$RfAtoAxis{z};
						my $tangle = rad2deg(acos($dotproduct/($AxAtoAxismag*$RfAtoAxismag)));
						unless ($tangle> 90) { #checking the test axis vector is in the same direction as the reference
							push @Axangles, $tangle;
							push @tAxisx, $AxAtoAxis{x};
							push @tAxisy, $AxAtoAxis{y}; 
							push @tAxisz, $AxAtoAxis{z};
						} else {
							push @notAxangles, $tangle;
						}
						my $tdev = (max @tAxisx - min @tAxisx) + (max @tAxisy - min @tAxisy) + (max @tAxisz - min @tAxisz);
						if (abs($tdev) < abs($dev)){
							$dev = $tdev;
							@Axisx = @tAxisx;
							@Axisy = @tAxisy;
							@Axisz = @tAxisz;
							$selaxi = $Axi;
						}
					}
				}
			}
		}
	}
}
my %AvAxis; #storing the average axis, this is an average of the normal vectors from Axi points (I hope)
unless (@Axisx and @Axisy and @Axisz){
	print <<ENDERROR;
Content-type: text/html

	
<!DOCTYPE HTML>
<html>
$head
<body><p>$warnings</p>
</body>
$foot
</html>
ENDERROR
die;
}
$AvAxis{x} = (sum @Axisx)/ (scalar @Axisx);
$AvAxis{y} = (sum @Axisy)/ (scalar @Axisy);
$AvAxis{z} = (sum @Axisz)/ (scalar @Axisz);
my $AvAxismag = sqrt(($AvAxis{x}**2)+($AvAxis{y}**2)+($AvAxis{z}**2));
my %Axisu; #unit vector for axis
$Axisu{x} = $AvAxis{x}/$AvAxismag;
$Axisu{y} = $AvAxis{y}/$AvAxismag;
$Axisu{z} = $AvAxis{z}/$AvAxismag;
my %CenAx; #average of Axi coordinates from chains
my %SumAxi;
for (my $i=0; $i < $mer; $i++) {
	$SumAxi{x} += $x{"$chains[$i]",$selaxi};
	$SumAxi{y} += $y{"$chains[$i]",$selaxi};
	$SumAxi{z} += $z{"$chains[$i]",$selaxi};
	}
$CenAx{x} = ($SumAxi{x} / $mer);
$CenAx{y} = ($SumAxi{y} / $mer);
$CenAx{z} = ($SumAxi{z} / $mer);
my $CenAxdotAxisu = ($CenAx{x}*$Axisu{x})+($CenAx{y}*$Axisu{y})+($CenAx{z}*$Axisu{z});
my @tilts;
my @swings;
my %HN;
my %HC;
my %HNtoHC;
my %HNtoAvAxis;
my $HNtoAvAxismag;
my $HNtoHCmag;
my $AxisdotHNtoHC;
foreach my $chain (@chains) { #iterate for each chain using same axis definition
	%HN; #calculate and store helix end locations as hashes (ends defined as average of first and last 4 residues)
	$HN{x} = ($x{$chain,$Ni}+$x{$chain,$Ni+1}+$x{$chain,$Ni+2}+$x{$chain,$Ni+3})/4;
	$HN{y} = ($y{$chain,$Ni}+$y{$chain,$Ni+1}+$y{$chain,$Ni+2}+$y{$chain,$Ni+3})/4;
	$HN{z} = ($z{$chain,$Ni}+$z{$chain,$Ni+1}+$z{$chain,$Ni+2}+$z{$chain,$Ni+3})/4;
	%HC;
	$HC{x} = ($x{$chain,$Ci}+$x{$chain,$Ci-1}+$x{$chain,$Ci-2}+$x{$chain,$Ci-3})/4;
	$HC{y} = ($y{$chain,$Ci}+$y{$chain,$Ci-1}+$y{$chain,$Ci-2}+$y{$chain,$Ci-3})/4;
	$HC{z} = ($z{$chain,$Ci}+$z{$chain,$Ci-1}+$z{$chain,$Ci-2}+$z{$chain,$Ci-3})/4;
	%HNtoHC; #calculate vector along helix to calculate tilts
	$HNtoHC{x} = $HC{x} - $HN{x};
	$HNtoHC{y} = $HC{y} - $HN{y};
	$HNtoHC{z} = $HC{z} - $HN{z};
	unless ($HNtoHC{x} + $HNtoHC{y} + $HNtoHC{z} != 0){
		print <<ENDERROR;
Content-type: text/html

	
<!DOCTYPE HTML>
<html>
$head
 <body>
 <p>
 $warnings
 </p>
 </body>
$foot
</html>
ENDERROR
 die;
		}
	$HNtoHCmag = sqrt(($HNtoHC{x}**2)+($HNtoHC{y}**2)+($HNtoHC{z}**2));
	$AxisdotHNtoHC = ($AvAxis{x}*$HNtoHC{x})+($AvAxis{y}*$HNtoHC{y})+($AvAxis{z}*$HNtoHC{z});
	my $tilt = rad2deg(acos($AxisdotHNtoHC/($AvAxismag*$HNtoHCmag)));
	push(@tilts, $tilt);
	my %CenAxtoHN;
	$CenAxtoHN{x} = $HN{x} - $CenAx{x};
	$CenAxtoHN{y} = $HN{y} - $CenAx{y};
	$CenAxtoHN{z} = $HN{z} - $CenAx{z};
	my %CenAxtoHC;
	$CenAxtoHC{x} = $HC{x} - $CenAx{x};
	$CenAxtoHC{y} = $HC{y} - $CenAx{y};
	$CenAxtoHC{z} = $HC{z} - $CenAx{z};
	my $CenAxtoHNdotAxisu = ($CenAxtoHN{x}*$Axisu{x})+($CenAxtoHN{y}*$Axisu{y})+($CenAxtoHN{z}*$Axisu{z});
	my $CenAxtoHCdotAxisu = ($CenAxtoHC{x}*$Axisu{x})+($CenAxtoHC{y}*$Axisu{y})+($CenAxtoHC{z}*$Axisu{z});
	my %CenN;
	$CenN{x} = $CenAx{x} + ($CenAxtoHNdotAxisu*$Axisu{x});
	$CenN{y} = $CenAx{y} + ($CenAxtoHNdotAxisu*$Axisu{y});
	$CenN{z} = $CenAx{z} + ($CenAxtoHNdotAxisu*$Axisu{z});
	my %CenC;
	$CenC{x} = $CenAx{x} + ($CenAxtoHCdotAxisu*$Axisu{x});
	$CenC{y} = $CenAx{y} + ($CenAxtoHCdotAxisu*$Axisu{y});
	$CenC{z} = $CenAx{z} + ($CenAxtoHCdotAxisu*$Axisu{z});
	my %NtoCen;
	$NtoCen{x} = $CenN{x} - $HN{x};
	$NtoCen{y} = $CenN{y} - $HN{y};
	$NtoCen{z} = $CenN{z} - $HN{z};
	my $NtoCenmag = sqrt(($NtoCen{x}**2)+($NtoCen{y}**2)+($NtoCen{z}**2));
	my %CentoC;
	$CentoC{x} = $HC{x} - $CenC{x};
	$CentoC{y} = $HC{y} - $CenC{y};
	$CentoC{z} = $HC{z} - $CenC{z};
	my %NtoCflat;
	$NtoCflat{x} = $CentoC{x} + $NtoCen{x};
	$NtoCflat{y} = $CentoC{y} + $NtoCen{y};
	$NtoCflat{z} = $CentoC{z} + $NtoCen{z};
	my $NtoCflatmag = sqrt(($NtoCflat{x}**2)+($NtoCflat{y}**2)+($NtoCflat{z}**2));
	my $NtoCflatdotNtoCen = ($NtoCflat{x} * $NtoCen{x}) + ($NtoCflat{y} * $NtoCen{y}) + ($NtoCflat{z} * $NtoCen{z});
	my $swingdet = ($NtoCflat{x}*$NtoCen{y}*$Axisu{z})+($NtoCflat{z}*$NtoCen{x}*$Axisu{y})+($NtoCflat{y}*$NtoCen{z}*$Axisu{x})-($NtoCflat{z}*$NtoCen{y}*$Axisu{x})-($NtoCflat{x}*$NtoCen{z}*$Axisu{y})-($NtoCflat{y}*$NtoCen{x}*$Axisu{z});
	my $swing = rad2deg(atan2($swingdet,$NtoCflatdotNtoCen));
	push (@swings, $swing);
	}
my $tiltav = (sum @tilts) / (scalar @tilts);
my $swingav = (sum @swings) / (scalar @swings);
foreach (@tilts) { #round to 2 dp
	$_ = sprintf("%.2f", $_);
	}
foreach (@swings) { #round to 2 dp
	$_ = sprintf("%.2f", $_);
	}
my $table = "<tr>\n<td>Chain</td>\n<td>Tilt Angle    </td>\n<td>Swing Angle</td>\n</tr>\n";
for (my $i=0; $i < $mer; $i++) {
	$table .= "<tr>\n<td>$chains[$i]</td>\n<td>$tilts[$i]</td>\n<td>$swings[$i]</td>\n</tr>\n";
	}
$tiltav = sprintf("%.2f", $tiltav);
$swingav = sprintf("%.2f", $swingav);
$table .= "<tr>\n<td>Average</td>\n<td>$tiltav</td>\n<td>$swingav</td>\n</tr>\n";
if ($warnings =~ /^<b>Warnings and Errors<\/b><br>$/) {
	undef $warnings;
	}
{
print <<ENDANS;
Content-type: text/html


<!DOCTYPE HTML>
<html>
$head
 <body>
 <p>$warnings</p>
 <p>Angles calculated for helix $Ni-$Ci in $filename using axis residue $selaxi</p>
 <table>
 $table
 </table>
</body>
$foot
</html>
ENDANS
}