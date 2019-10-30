package  SoftwareConfig;
use Exporter 'import';
our @EXPORT = qw(parse_config);
##parse the software_files.config , and check the existence of each software or file
################################################

sub parse_config{
        my ($config,$soft)=@_;
        open IN,$config || die "Cannot open $config";
        my %ha;
        while(<IN>){
                chomp;
                next if /^#|^$/;
                s/\s+//g;
                my @p=split/=/,$_;
                $ha{$p[0]}=$p[1];
        }
        close IN;
        if(exists $ha{$soft}){
                if(-e $ha{$soft}){
                        return $ha{$soft};
                }else{
                         die "\nConfig Error: $soft path in $config does not exist\n";
                }
        }else{
                die "\nConfig Error: No default path to $soft in $config\n";
        }
}
1;              
__END__         
