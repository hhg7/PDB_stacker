#!/usr/bin/env perl

use strict;
use warnings;
use Math::Trig;#necessary for arcsin and arccos
use List::Util qw(min);#evaluates more quickly than a traditional if bracket, necessary for larger files

local $SIG{__WARN__} = sub {#kill the program if there are any warnings
	my $message = shift;
	print "$message\n";
	die;
};

sub Magnitude {
   my $fr = 0;
   foreach my $e (@{$_[0]}) {
      $fr += $e**2;
   }
   return($fr**0.5);
}

sub Cross_Product_3D {#takes the first 3 values of two arrays and returns 3 values
   my @fu = @{$_[0]};
   my @fv = @{$_[1]};
   my @r = ($fu[1]*$fv[2]-$fu[2]*$fv[1],-($fu[0]*$fv[2]-$fu[2]*$fv[0]),$fu[0]*$fv[1]-$fu[1]*$fv[0]);
   return($r[0],$r[1],$r[2]);
}

use constant RAD => atan2(1.0,1.0)/45.0;

my $OMEGA_MIN = 25;
my $DISTANCE_MIN = 4;
my $DISTANCE_MAX = 5.0;#ab initio limit for benzene stacking (J. Am. Chem. Soc. 124: 10887)
my $OMEGA_MAX = 50;
#my $DISTANCE_SLOPE = 1/($DISTANCE_MIN-$DISTANCE_MAX);
#my $DISTANCE_B = 1-$DISTANCE_SLOPE*$DISTANCE_MIN;
my $OMEGA_SLOPE = 1/($OMEGA_MIN-$OMEGA_MAX);#for y = mx + b
my $OMEGA_B = 1-$OMEGA_SLOPE*$OMEGA_MIN;#for y = mx + b

sub stacking_score {#this is set so I can easily change the values later. This takes three values: d, Xi, and Omega in that order.
   my $stack_score = 0;
   if ($_[0] <= $DISTANCE_MAX) {#two bases are not stacked if they are > 4.96 Angstroms away from one another
      if ($_[0] <= $DISTANCE_MIN) {
         $stack_score++;
#      } elsif (($_[0] >= $DISTANCE_MIN) && ($_[0] < $DISTANCE_MAX)) {
#         $stack_score += $DISTANCE_SLOPE*$_[0]+$DISTANCE_B;#for y = mx + b
      } elsif ($_[0] > $DISTANCE_MIN) {
         $stack_score += 1/($_[0]-$DISTANCE_MIN+1)**3;
#         print "$_[0] $stack_score\n";
      }
      if ($_[2] <= $OMEGA_MIN) {
         $stack_score++;
      } elsif (($_[2] > $OMEGA_MIN) && ($_[2] <= $OMEGA_MAX)) {
         $stack_score += $OMEGA_SLOPE*$_[2]+$OMEGA_B;#for y = mx + b
      } elsif ($_[2] > $OMEGA_MAX) {
         $stack_score = 0;
      }
      if (($_[1] > 45) && ($_[1] <= 135)) {#Xi only changes the the stacking score's sign
         $stack_score *= -1;#a T stack
      }
   }
   $stack_score = 100*$stack_score/2;#to get it into a %
   if ($stack_score > 100) {
      print "Entering d = $_[0], Xi = $_[1], and Omega = $_[2] gives $stack_score > 100\n";
      die;
   }
   return($stack_score);
}

print "This program examines AMBER output PDB files and outputs stacking data into files C1.A2.stack, etc.\n";

my (%results,%nucs,%MASS,%v);

my %atoms = (
   'A' => ['N1','C2','N3','C4','C5','C6','N6','N7','C8','N9'],
   'C' => ['N1','C2','O2','N3','C4','N4','C5','C6'],
   'G' => ['N1','C2','N3','C4','C5','C6','O6','N7','C8','N9'],
   'U' => ['N1','C2','O2','N3','C4','O4','C5','C6']
);

$MASS{A} = 2*(7+6+7+6+6+6+7+7+6+7);$MASS{C} = 2*(7+6+8+7+6+7+6+6);
$MASS{G} = 2*(7+6+7+6+6+6+8+7+6+7);$MASS{U} = 2*(7+6+8+7+6+8+6+6);

opendir(DIR,"pdbs") or die $!;#this is a 3x faster than glob
while (my $pdb = readdir(DIR)) {
	if ($pdb =~ m/^\.+$/) {
		next;
	}
   my $time;
   if ($pdb =~ m/(\d+)\.(\d)ns\.pdb$/) {
      $time = "$1.$2";
   } elsif ($pdb =~ m/(\d+)ns\.pdb$/) {
      $time = $1;
   } else {
      print "$pdb is not a pdb file.\n";
      next;
   }
   open(FH,"<pdbs/$pdb") or die "cannot read pdbs/$pdb: $!";
   $time = sprintf("%.1f",$time);
   if (!defined $time) {
       print "$pdb has no time.\n";
       die;
   }
  	my (%x,%y,%z,%magnitude); 
#      if (/^ATOM\s+\d+\s+([CNO])(\d)('?)\s+R?([ACGU])[35]?\s+(\d)\s+(-?)(\d+)\.(\d+)\s+(-?)(\d+)\.(\d+)\s+(-?)(\d+)\.(\d+)/) {#specific for AMBER PDB files, ignoring H atoms
   while (<FH>) {
   	if (/^ATOM\s+\d+\s+([CNO])(\d)('?)\s+R?([ACGU])[35]?\s+\S?\s+(\d+)\s+(-?)(\d+)\.(\d+)\s+(-?)(\d+)\.(\d+)\s+(-?)(\d+)\.(\d+)/) {
          $x{"$5 $1$2$3"} = "$6$7.$8";#nucleotide #, atom name, x coordinate, e.g. "3 H1'"
          $y{"$5 $1$2$3"} = "$9$10.$11";#nucleotide #, atom name, y coordinate
          $z{"$5 $1$2$3"} = "$12$13.$14";#nucleotide #, atom name, z coordinate
          $nucs{$5} = $4;#Residue map will be the same for every model
      }
   }
   close FH;
   my (%Center_of_mass,%normal,@u,@v);#normal = orthogonal
   foreach my $nuc (keys(%nucs)) {#Calculate center of mass and orthogonal=normal vectors for each nucleotide. More efficient than before
      $Center_of_mass{$nuc}{x} = 0;      $Center_of_mass{$nuc}{y} = 0;      $Center_of_mass{$nuc}{z} = 0;
      foreach my $atom (@{ $atoms{$nucs{$nuc}} }) {
         my $weight;
         if ($atom =~ m/C/) {
            $weight = 12;
         } elsif ($atom =~ m/N/) {
            $weight = 14;
         } elsif ($atom =~ m/O/) {
            $weight = 16;
         } else {
            print "Could not get a weight.\n";
            die;
         }
         $Center_of_mass{$nuc}{x} += $weight*$x{"$nuc $atom"};
         $Center_of_mass{$nuc}{y} += $weight*$y{"$nuc $atom"};
         $Center_of_mass{$nuc}{z} += $weight*$z{"$nuc $atom"};
      }
      $Center_of_mass{$nuc}{x} /= $MASS{$nucs{$nuc}};
      $Center_of_mass{$nuc}{y} /= $MASS{$nucs{$nuc}};
      $Center_of_mass{$nuc}{z} /= $MASS{$nucs{$nuc}};
      if ($nucs{$nuc} eq 'C') {#https://en.wikipedia.org/wiki/Cross_product,0s are place holders to match their computation notation
#1st intrabase vectors should have normal vector ABOVE the plane of the base
#new origin is from center of mass of current nucleotide for EVERYTHING.
         @u = ($x{"$nuc O2"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc O2"}-$Center_of_mass{$nuc}{y},$z{"$nuc O2"}-$Center_of_mass{$nuc}{z});
         @v = ($x{"$nuc N4"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc N4"}-$Center_of_mass{$nuc}{y},$z{"$nuc N4"}-$Center_of_mass{$nuc}{z});
      } elsif ($nucs{$nuc} eq 'U') {
         @u = ($x{"$nuc O2"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc O2"}-$Center_of_mass{$nuc}{y},$z{"$nuc O2"}-$Center_of_mass{$nuc}{z});  
         @v = ($x{"$nuc O4"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc O4"}-$Center_of_mass{$nuc}{y},$z{"$nuc O4"}-$Center_of_mass{$nuc}{z});
      } elsif ($nucs{$nuc} eq 'A') {
         @u = ($x{"$nuc N6"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc N6"}-$Center_of_mass{$nuc}{y},$z{"$nuc N6"}-$Center_of_mass{$nuc}{z});
         @v = ($x{"$nuc C8"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc C8"}-$Center_of_mass{$nuc}{y},$z{"$nuc C8"}-$Center_of_mass{$nuc}{z});
      } elsif ($nucs{$nuc} eq 'G') {
         @u = ($x{"$nuc O6"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc O6"}-$Center_of_mass{$nuc}{y},$z{"$nuc O6"}-$Center_of_mass{$nuc}{z});
         @v = ($x{"$nuc C8"}-$Center_of_mass{$nuc}{'x'},$y{"$nuc C8"}-$Center_of_mass{$nuc}{y},$z{"$nuc C8"}-$Center_of_mass{$nuc}{z});
      }
      ($normal{$nuc}{above}{x},$normal{$nuc}{above}{y},$normal{$nuc}{above}{z}) = &Cross_Product_3D(\@u,\@v);
      ($normal{$nuc}{below}{x},$normal{$nuc}{below}{y},$normal{$nuc}{below}{z}) = &Cross_Product_3D(\@v,\@u);
#Calculate lengths of normal vectors.  Magnitude above and below are equal, so there is no point in caculating both.
      @u = ($normal{$nuc}{below}{x},$normal{$nuc}{below}{y},$normal{$nuc}{below}{z});
      $magnitude{$nuc} = &Magnitude(\@u);
      foreach my $e (keys($normal{$nuc})) {#above and below.  Need to add coordinates back after coordinates removed during @u and @v declaration above
         foreach my $f (keys($normal{$nuc}{$e})) {#x, y, and z
            $normal{$nuc}{$e}{$f} += $Center_of_mass{$nuc}{$f};#6 lines are now 4 through these two foreach loops.
         }
      }
   }#done getting normal vectors and centers of mass for each nucleotide
   foreach my $nuc0 (keys(%nucs)) {
      foreach my $nuc1 (keys(%nucs)) {
         if ($nuc0 >= $nuc1) {#only calculate 1-2, not 2-1, etc.
            next;
         }
         my $d0 = (($Center_of_mass{$nuc0}{x}-$Center_of_mass{$nuc1}{x})**2#distance between centers of mass
   		+($Center_of_mass{$nuc0}{y}-$Center_of_mass{$nuc1}{y})**2
		+($Center_of_mass{$nuc0}{z}-$Center_of_mass{$nuc1}{z})**2)**0.5;
         $v{d}{$nuc0}{$nuc1}{$time} = $d0;
#         if ($d0 > $DISTANCE_MAX) {#don't waste time calculating omega and Xi if they're discounted
#            $results{$nuc0}{$nuc1}{$time} = 0;
#            next;
#         }
         my $d2 = min((($normal{$nuc0}{below}{x}-$Center_of_mass{$nuc1}{x})**2+($normal{$nuc0}{below}{y}-$Center_of_mass{$nuc1}{y})**2+($normal{$nuc0}{below}{z}-$Center_of_mass{$nuc1}{z})**2)**0.5,(($normal{$nuc0}{above}{x}-$Center_of_mass{$nuc1}{x})**2+($normal{$nuc0}{above}{y}-$Center_of_mass{$nuc1}{y})**2+($normal{$nuc0}{above}{z}-$Center_of_mass{$nuc1}{z})**2)**0.5);#is the vector above the plane or below the plane closer to the next nucleotide's center of mass?
         my $omega = (1/RAD)*acos( ($d0*$d0+$magnitude{$nuc0}*$magnitude{$nuc0}-$d2*$d2)/(2*$d0*$magnitude{$nuc0}) );
         $v{omega}{$nuc0}{$nuc1}{$time} = $omega;
         my $Xi = 387420489;#impossible of course, just as something too big
         @u = ($normal{$nuc0}{above}{x}-$Center_of_mass{$nuc0}{x},$normal{$nuc0}{above}{y}-$Center_of_mass{$nuc0}{y},$normal{$nuc0}{above}{z}-$Center_of_mass{$nuc0}{z});
         foreach my $v ('above','below') {#pretend that both vectors are centered on the 5' nucleotide's Center of Mass
            @v = ($normal{$nuc1}{$v}{x}-$Center_of_mass{$nuc0}{x},$normal{$nuc1}{$v}{y}-$Center_of_mass{$nuc0}{y},$normal{$nuc1}{$v}{z}-$Center_of_mass{$nuc0}{z});
#            print "<$u[0],$u[1],$u[2]> x <$v[0],$v[1],$v[2]>\n";
            my @tmpCrossProduct = &Cross_Product_3D(\@u,\@v);
            my $tmpXi = (1/RAD)*asin((&Magnitude(\@tmpCrossProduct))/(&Magnitude(\@u)*&Magnitude(\@v)));
            if ($tmpXi < $Xi) {
               $Xi = $tmpXi;
            }
         }
         $v{Xi}{$nuc0}{$nuc1}{$time} = $Xi;
         if ($omega > $OMEGA_MAX) {#don't waste time score Xi if it's discounted anyway
            $results{$nuc0}{$nuc1}{$time} = 0;
            next;
         }
         $results{$nuc0}{$nuc1}{$time} = sprintf("%.1f",&stacking_score($d0,$Xi,$omega));
#         if ((49.7 < $results{$nuc0}{$nuc1}{$time}) && ($results{$nuc0}{$nuc1}{$time} < 50.5)) {
#            print "$omega $d0 => $results{$nuc0}{$nuc1}{$time}\n";
#         }
#         print "$nuc0-$nuc1 $time ns = $results{$nuc0}{$nuc1}{$time}\n";
      }#for each $nuc1
   }#for each $nuc0
}#end of while loop
closedir(DIR);
print "Finished reading PDBs\n";
foreach my $nuc0 (keys(%results)) {
   foreach my $nuc1 (keys($results{$nuc0})) {
      if ($nuc0 >= $nuc1) {
         next;
      }
      my $name = "$nucs{$nuc0}$nuc0-$nucs{$nuc1}$nuc1" . "_o$OMEGA_MIN-$OMEGA_MAX" . "_d$DISTANCE_MIN" . "_cube_stack";
      open(STACK,">$name") or die "cannot write to $name: $!";
      open(OMEGA,">OMEGA_$nucs{$nuc0}$nuc0-$nucs{$nuc1}$nuc1-stack") or die $!;
      open(FD,">D_$nucs{$nuc0}$nuc0-$nucs{$nuc1}$nuc1-stack") or die $!;
      open(XI,">Xi_$nucs{$nuc0}$nuc0-$nucs{$nuc1}$nuc1-stack") or die $!;
      foreach my $time (sort {$a <=> $b} keys($results{$nuc0}{$nuc1})) {#for each timepoint
         print STACK "$time $results{$nuc0}{$nuc1}{$time}\n";
         print OMEGA "$time $v{omega}{$nuc0}{$nuc1}{$time}\n";
         print FD "$time $v{d}{$nuc0}{$nuc1}{$time}\n";
         print XI "$time $v{Xi}{$nuc0}{$nuc1}{$time}\n";
      }
      close STACK;     close OMEGA;     close FD;      close XI;
      print "Wrote files $name, OMEGA.$nucs{$nuc0}$nuc0.$nucs{$nuc1}$nuc1-stack, D.$nucs{$nuc0}$nuc0.$nucs{$nuc1}$nuc1-stack, and Xi.$nucs{$nuc0}$nuc0.$nucs{$nuc1}$nuc1-stack\n";
   }
}
