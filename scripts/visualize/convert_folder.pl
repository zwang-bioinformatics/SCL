my $dir = $ARGV[0]."/";
my $dir2 = $ARGV[1]."/";
my $length = $ARGV[2];
my $convert = $ARGV[3];
opendir(DIR,$dir);
my @files = readdir(DIR);
my $count=1;
my $dir1=$ARGV[1];
`mkdir $dir1/`;
foreach my $file (@files){
        if($file eq '.' || $file eq '..'|| $file eq " " || $file eq''){
                next
        }
	my @item = split(/\./,$file);
	$name = $dir2.$item[0].".pdb";
        #$name = $dir.$count.".pdb";
	#print $name."\n";
        $name1 = $dir.$file;
        `perl convert_single.pl $name1 $name $length $convert`;
	#`cp $name $dir1/$name`;
		
        $count++;
}
#for(my $i=1;$i<$count;$i++){
#	my $name = $dir.$i.".pdb";
#	my $name2 = "$dir1/".$i.".pdb"."\n";
#	`cp $name $name2`;
#	`rm $name`;
#}
