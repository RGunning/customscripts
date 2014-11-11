#!/usr/bin/env perl
####################################################
# Created By: Nic Walker (nw11@sanger.ac.uk)       #
# Modified By: Richard Gunning (rg12@sanger.ac.uk  #
#												   #
####################################################


# getopt
# json
# read the vocab file
# get out the terms for gender into a hash
# build up a hash for each track
# track->{sex}, track->{assay}, track->{}

# test:
# perl bin/make-wueb-metadata --vocab_url http://jeremy-bio.sysbiol.cam.ac.uk/EPIGATEWAY/blueprint/bp-vocab.json --url_base 'http://a.b.c.com' --url_suffix '.gz'  --sexes 'M F M' --assays 'OX BS BS' --cell_types 'B T B' --urls 'thing1.gz thing2.gz thing3.gz'
# perl bin/make-wueb-metadata --vocab_url http://jeremy-bio.sysbiol.cam.ac.uk/EPIGATEWAY/blueprint/bp-vocab.json  --sexes 'M F M' --assays 'OX BS BS' --cell_types 'B T B' --url_names 'thing1.gz thing2.gz thing3.gz' --base_url 'http://a.b.c.com' --url_suffix '.gz'  > out1

use Modern::Perl;
use Getopt::Long;
use Data::Dumper;
use Path::Tiny;
use JSON;
use HTTP::Tiny;
use Storable 'dclone';
use List::AllUtils qw(uniq);
use Set::CrossProduct;
use Color::Scheme;

print STDERR Dumper \@ARGV;

my( $file, $sexes, $assays,$url_names, $base_url, $help, $json_url, $cell_types, $url_suffix, $hmark, $strain, $repeat, $peak);
my $name_suffix;
my $colour_by = 'cell_types';
my $type = 'bedgraph';
GetOptions( "vocab_url=s"   => \$json_url,
			"strain=s"		=> \$strain,
            "sexes=s"    	=> \$sexes,
            "cell_types=s" 	=> \$cell_types,
            "assays=s"   	=> \$assays,
            "hmark=s"		=>	\$hmark,
            "url_names=s"   => \$url_names,
            "url_suffix=s" 	=> \$url_suffix,
            "name_suffix=s" => \$name_suffix,
            "base_url=s" 	=> \$base_url,
            "colour_by=s" 	=> \$colour_by,
            "type=s"        => \$type,
            "repeat=s"		=> \$repeat,
            "peak=s"		=> \$peak,
            "help" => \$help
            )
    or die("Error in command line arguments\n");



if( $help ){

my $help_message = <<HELP;
   make-wueb-metadata
   Usage:
     make-wueb-metadata --vocab_url --base_url --sexes --cell_types --assays

   E.G.
   make-wueb-metadata --vocab_url http://jeremy-bio.sysbiol.cam.ac.uk/EPIGATEWAY/blueprint/bp-vocab.json  --base_url 'http://a.b.c.com' --url_suffix '.gz' --sexes 'M F M' --assays 'OX BS BS' --cell_types 'B T B' --url_names 'thing1.gz thing2.gz thing3.gz' --type bedgrpah
HELP

print $help_message;
exit 0;
}

die "url_names not defined"       if ! defined($url_names);
die "url_suffix not defined"       if ! defined($url_suffix);
die "url_base not defined"   if ! defined($base_url);


my @colour_by = split /\s/, $colour_by;

#my @inputs = ( $strain, $sexes, $assays, $cell_types, $hmark, $repeat);

# translate strings into arrays
#my @properties;

my @url_names = map { $_ . $url_suffix } split /\s/, $url_names;
#push @properties, \@url_names;
#foreach my $input (@inputs){
#   my @input = split /\s/, $input;
#   push @properties, \@input;
#}

print STDERR "Setting Properties:\n";

# SET PROPERTY HASH BASED ON WHAT USER SUPLIES
my %properties = (
    url_names  	=>  \@url_names,
);

$properties{'strain'}= [split /\s/, $strain] if defined($strain);
$properties{'sexes'} = [split /\s/, $sexes] if defined($sexes);
$properties{'assays'} = [split /\s/, $assays] if defined($assays);
$properties{'cell_types'} = [split /\s/, $cell_types] if defined($cell_types);
$properties{'hmark'} = [split /\s/, $hmark] if defined($hmark);
$properties{'repeat'} = [split /\s/, $repeat] if defined($repeat);
$properties{'peak'} = [split /\s/, $peak] if defined($peak);

print STDERR Dumper \%properties;
print STDERR "====\n";

# download vocab_url and make vocab hash
# called term_id_hash, has e.g. { OX => '3', BS=>'4'  'T' =>'2' 'B' => 1}
my $content = HTTP::Tiny->new->get( $json_url );
die "Failed!\n" unless $content->{success};
$content=$content->{content};
print STDERR "Got vocab json file:\n";
print STDERR Dumper $content;
print STDERR "====\n";

my $vocab_hash = from_json( $content );
my $terms = $vocab_hash->{terms};
#print Dumper $terms;
my %term_id_hash;
foreach my $term_key ( keys %$terms ) {
        my $term_array = $terms->{$term_key};
        my $term_id = $term_array->[0];
        $term_id_hash{$term_id} = $term_key;
        #print STDERR "$term_key: TERMID: $term_id\n";
}

print STDERR "TERM_ID_HASH\n";
print STDERR Dumper \%term_id_hash;
print STDERR "====\n";

my @datahub;

my $track_template;

if( $type eq 'bedgraph'){
    $track_template = {
      type => "bedgraph",
      mode => "show",
      colorpositive => "#ff33cc/#B30086",
      backgroundcolor => "#ffe5ff",
      metadata => { blueprint => [] },
      height => 40,
      group  => 1,
      name   => ''
  };
}

if( $type eq 'hammock' ){
  $track_template = {
      type => "hammock",
      mode => "barplot",
      metadata => { blueprint => [] },
      showscoreidx => "0",
	  scorenamelst => ["signal value", "P value (-log10)","Q value (-log10)"]
  };
}

my $property_iter = 0;
my $item_iter = 0;

say STDERR "COLOUR SCHEME";
my $properties2 = dclone \%properties;
delete $properties2->{url_names};
   # sort on the tmp col
  # my @c = sort { $a->[ $sort_col ] cmp $b->[ $sort_col] } @b;
  # # remove the tmp sort col
  # my @d = map { [@{ $_ }[ @all_cols ] ] } @c;
  # return \@d;
my $colour_scheme = color_scheme_factory( $properties2, @colour_by );


say STDERR Dumper $colour_scheme;
print STDERR "===========\n" ;



#foreach my $property ( @properties){
# properties are cell_type, sex, assay, url_name ...
# i.e. these are columns in the datasource
foreach my $property_key ( sort keys %properties ){

   my $property = $properties{$property_key};
   say STDERR "Processing: $property_key";
   ## This is is essentially going through value of the property
   ## e.g. for sex: M F M F
   ## which correponds to the row in the datasource column
   $item_iter = 0;
   foreach my $property_value ( @$property ){
       #print "assign $property_value in $item_iter\n";

     if($property_iter == 0){
          my $track = dclone $track_template ;
          $datahub[$item_iter] = $track;
     }

     if( $property_key eq "url_names" ) {
         # url property, and assign the template
         #my $track = dclone $track_template ;
         #$track->{url} = $base_url . "/" . $property_value ;
         $datahub[$item_iter]->{url} = $base_url . "/" . $property_value ; #$track;
         if ($property_value =~ /broad/){
         	$datahub[$item_iter]->{scorenamelst}=["signal value", "P value (-log10)"];
         }
     } else {
         #get the correct reference
         my $term_num = $term_id_hash{$property_value};
         my $metadata = $datahub[$item_iter]->{metadata}->{blueprint};
         push @$metadata, $term_num;
         my @name;
         my $prev_name = $datahub[$item_iter]->{name};
         say STDERR "previous name ", $prev_name;
         if( $prev_name =~ /-/ ){
           @name = split /-/, $prev_name ;
         } elsif ($prev_name ne '' ){
           push @name, $prev_name;
         }

         say STDERR "GOT the previous name in array: @name";
         push @name, $property_value;
         say STDERR "pushed $property_value onto name array";
         my $name;
         if( @name > 0 ){
             $name = join "-", @name;
         }else{
             $name = $name[0];
         }

         say STDERR "setting $name to name";
         $datahub[$item_iter]->{name} = $name;
         $datahub[$item_iter]->{colorpositive} = $colour_scheme->{ $name };
     }

     if ($property_iter == scalar( keys %properties) -1 ){
         if (defined( $name_suffix ) ){
             $datahub[$item_iter]->{name} .= $name_suffix;
         }
     }

     $item_iter++;



   }
   $property_iter++;
}

## put on metadata references
my $metadata_entry = {
   type =>  "metadata",
   vocabulary_set =>  {
            "blueprint" => "$json_url",
    }
};

push @datahub, $metadata_entry;

print STDERR "DATAHUB:\n";
print STDERR Dumper @datahub;
print STDERR "====\n";

my $json = to_json(\@datahub, {pretty=> 1});
print STDERR "JSON PRINTED TO STDOUT\n";
print $json;

sub color_scheme_factory {

   my $properties_hash = shift;
   my %properties_to_colour = map { $_ => 1 } @_;

   say STDERR "properties_hash";
   say STDERR Dumper $properties_hash;
   say STDERR Dumper %properties_to_colour;

   # get a list of colours
   my $colours = get_colours();
   print STDERR "list of colours:\n";
   print STDERR Dumper $colours;

   # filter the properties hash
   #my @keys = grep  { $properties_to_colour{$_} } sort keys %$properties_hash;
   #say STDERR "only need @keys";
   #die "No valid property for colouring" if @keys == 0;
   # get the groups of differnt unique values
   my @values;
   my $iter = 0;
   my @sort_cols;
   foreach my $property (sort keys %$properties_hash){
       my $vals = $properties_hash->{ $property };
       my @unique_vals = uniq @$vals;
       push @values, \@unique_vals;
       if( exists $properties_to_colour{ $property } ){
           push @sort_cols, $iter;
       }
       $iter++;
   }
   my $s = Set::CrossProduct->new(\@values);
   my $combinations = $s->combinations;

   # sort by the columns
   my $colour_combinations = colour_AoA($combinations, @sort_cols );
   say STDERR "Coloured Combinations";
   say STDERR Dumper $colour_combinations;

   return $colour_combinations;
}

sub get_colours{

   my $scheme = Color::Scheme->new
      ->from_hex('B23737')
      ->scheme('analogic')
      ->distance(0.3)
      ->add_complement(1)
      ->web_safe(1);

   my @list = $scheme->colors();
   return \@list;
}


sub colour_AoA {
   my $AoA = shift;
   my @sort_cols = @_;
   my $num_cols = scalar( @{  $AoA->[0] } );
   my $sort_col = $num_cols; # since start from 0
   my @all_cols = ( 0 .. $num_cols - 1 );

   # make temp sort col
   my @b = map {  [ @{$_}, join '-',  @{ $_ }[ @sort_cols ] ] } @$AoA;

   #make a colour hash
   my @T = map { join '-', @{ $_ }[ @sort_cols ]  }  @$AoA;
   my $colour_list = get_colours();
   my %colour_hash;
   my $iter = 0;
   foreach my $t (@T) {
       if( ! exists($colour_hash{ $t } )){
           $colour_hash{$t} = '#' . $colour_list->[$iter];
           $iter++;
       }
   }
   say STDERR "COLOUR_HASH IN colour_AoA";
   say STDERR Dumper %colour_hash;
   my %colour_map = map {
                        join( '-', @{$_}[ @all_cols ])
                         =>
                         $colour_hash{ $_->[$sort_col] }

                         } @b;

  return \%colour_map;
}



#print STDERR "COLOUR SCHEME\n";
#print STDERR Dumper color_scheme_factory( \%properties, 'cell_types' );
#  submit values
#  gives back a hash
#  that allows colouring by multiple properties
#  e.g. M-OX-BS if for three properties. This can
#  be used above to check the name
# properties => (value,value,value) is e.g. sex => (M, F, M,)
# =cut
# sub color_scheme_factory_old {
#
#    my $properties_hash = shift;
#    my %properties_to_colour = map { $_ => 1 } @_;
#
#    print STDERR "properties_hash";
#    print STDERR Dumper $properties_hash;
#    print STDERR Dumper %properties_to_colour;
#
#    # get a list of colours
#    my $colours = get_colours();
#    print STDERR "list of colours:\n";
#    print STDERR Dumper $colours;
#
#    # filter the properties hash
#    #my @keys = grep  { $properties_to_colour{$_} } sort keys %$properties_hash;
#    #say STDERR "only need @keys";
#    #die "No valid property for colouring" if @keys == 0;
#    # get the groups of differnt unique values
#    my @values;
#    foreach my $property (sort keys %$properties_hash;){
#        my $vals = $properties_hash->{ $property };
#        my @unique_vals = uniq @$vals;
#        push @values, \@unique_vals;
#    }
#
#    say STDERR "uniq values";
#    say STDERR Dumper @values;
#
#    # get the different combinations if more than one property
#    my %colour_hash;
#    if( @values > 1){
#        my $s = Set::CrossProduct->new(\@values);
#        my $combinations = $s->combinations;
#        # now join these into a hash and assign a color
#        foreach my $combination (@$combinations){
#            foreach my $el (@$combination){
#                $el
#            }
#             my $property_str = join "-", @$combination;
#             $colour_hash{$property_str} = shift @$colours;
#        }
#     } else{
#        my $first_values = shift @values;
#        foreach my $value (@$first_values){
#           $colour_hash{$value} = shift @$colours;
#        }
#     }
#     return \%colour_hash;
# }
#
# # FIX THIS
# # THIS DOES NOT WORK, YOU HAVE TO ENUMERATE THE WHOLE HASH AND COLOUR EACH COMBINATION
# # SPECIFIED BY CELL TYPE OR WHATEVER.
#
# sub assign_colour_hash {
#      my $AoA = shift;
#      my $num_sort_cols = shift;
#      my $total_cols = scalar ( @{  $AoA } );
#      say STDERR "In assign colour hash";
#      say STDERR "Number of sort cols = $num_sort_cols";
#      say STDERR "Number of total cols = $total_cols";
#      # get colours
#      my %combination_colour_hash;
#      my $colours = get_colours();
#      my $idx = 0;
#      my $iter;
#      foreach my $array (@$AoA){
#          my $str = join "-", @$array;
#          $combination_colour_hash{$str} = $colours->[$idx];
#          my $idx = int( $iter /  ( $total_cols / $num_sort_cols ) );
#      }
#      return \%combination_colour_hash;
# }sub assign_colour_hash {
#      my $AoA = shift;
#      my $num_sort_cols = shift;
#      my $total_cols = scalar ( @{  $AoA } );
#      say STDERR "In assign colour hash";
#      say STDERR "Number of sort cols = $num_sort_cols";
#      say STDERR "Number of total cols = $total_cols";
#      # get colours
#      my %combination_colour_hash;
#      my $colours = get_colours();
#      my $idx = 0;
#      my $iter;
#      foreach my $array (@$AoA){
#          my $str = join "-", @$array;
#          $combination_colour_hash{$str} = $colours->[$idx];
#          my $idx = int( $iter /  ( $total_cols / $num_sort_cols ) );
#      }
#      return \%combination_colour_hash;
# }
#
#
# sub sort_AoA {
#    my $AoA = shift;
#    my @sort_cols = @_;
#    my $num_cols = scalar( @{  $AoA->[0] } );
#    my $sort_col = $num_cols; # since start from 0
#    my @all_cols = ( 0 .. $num_cols - 1 );
#
#    # make temp sort col
#    my @b = map {  [ @{$_}, join '-',  @{ $_ }[ @sort_cols ] ] } @$AoA;
#    # sort on the tmp col
#    my @c = sort { $a->[ $sort_col ] cmp $b->[ $sort_col] } @b;
#    # remove the tmp sort col
#    my @d = map { [@{ $_ }[ @all_cols ] ] } @c;
#    return \@d;
# }
# =cut
#

