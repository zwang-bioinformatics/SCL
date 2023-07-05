my $dir= $ARGV[0];
my $file1=$ARGV[1];
#$my $dir1=$ARGV[1];
opendir(DIR,$dir);
open(OUT,">$file1");
my @files = readdir(DIR);
my $count=1;
#`mkdir $dir1`;
my @name_f;
foreach my $file (@files){
	if($file eq '.' || $file eq '..'|| $file eq " " || $file eq''){
		next
	}
	#print  $file."\n";		
	$name_f[$count-1]=$file;
	$count++;
}
#print $count."\n";
my %hash;
my @array;
my @array1;
my @a;
for(my $i=1;$i<$count;$i++){
	my @a2;
	for(my $j=1;$j<$count;$j++){
		$a2[$j-1]=0;
		#print $a2[$j-1];	
	}
	push @a,[@a2];
	#print "\n";
}

for(my $i=1;$i<$count;$i++){
	my $sum = 0;
	for(my $j=1;$j<$count;$j++){
		if($i < $j){
			my $name = $dir.$name_f[$i-1];
			my $name2 = $dir.$name_f[$j-1];
			my $name3 = $i."_$j.txt";
			my $name4 = $i."_$j.txt_atm";
			`./TMscore $name $name2 -o $name3`;
			`rm $name3`;
			open(IN,"$name4");
			while(<IN>){
				my $line = $_;
				if(substr($line,8,1) eq "T"){
					my $score = substr($line,17,6);
					$a[$i-1][$j-1]=$score;
					$sum = $sum + $score;	
				}
			}
			close(IN);
			`rm $name4`;
			#print $j;
		}
		if($i>$j){
			my $score = $a[$j-1][$i-1];
			$sum = $sum + $score;
			#print $j;
		}
	}
	
	$array[$i]=$sum/($count-1);
	#print "\n".$i."th qscore is ".$array[$i]."\n";
	$array1[$i]=$i;
	$hash{$i}=$array[$i];
}

my $num=1;
#for(my $i=1;$i<$count;$i++){
	#my $name = $dir.$i.".pdb";
	#my $name2 = "$dir1/".$i.".pdb"."\n";
	#`cp $name $name2`;
	#`rm $name`;
#}

foreach my $name (sort{$hash{$b} <=> $hash{$a}} keys %hash){
	if($num<=3){
		print "top ".$num." ".$name_f[$name-1]."	".$hash{$name}."\n";
	}
	$num++;
	print OUT $name_f[$name-1]." qscore: ".$hash{$name}."\n";
}
close(OUT);
