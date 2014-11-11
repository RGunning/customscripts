#!/usr/bin/env perl
#===============================================================================
#
#         FILE: lsf_job_manager.pl
#
#        USAGE: ./lsf_job_manager.pl
#
#  DESCRIPTION:
#
#      OPTIONS: ---
# REQUIREMENTS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: NJWALKER (), nw11@sanger.ac.uk modified RJGUNNING, rg12@sanger.ac.uk
#      COMPANY:
#      VERSION: 1.0
#      CREATED: 24/02/13 12:05:44
#     REVISION: ---
#===============================================================================

use strict;
use warnings;
use strict;
use warnings;
use YAML::Any;
use Graph::Easy;
use Graph;
use Graph::Convert;
use Path::Class qw(dir file);
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::Job;
use Time::Stamp qw(-stamps);
use autobox::Core;
# This should just be given the directory as an input
# Then everything else should be down to reading the config file that
# gives the names of the appropriate files.

my $pipeline_dir = dir($ARGV[0]);
print "$pipeline_dir\n";
my $pipeline_name = $pipeline_dir->basename;
print "$pipeline_name\n";
$pipeline_name =~ s/.graph.yaml//g;
my $pipeline_graph_and_jobs_yaml_file = file($pipeline_dir);
print "$pipeline_graph_and_jobs_yaml_file\n";
my $pipeline_graph_and_jobs_yaml = $pipeline_graph_and_jobs_yaml_file->slurp;
my $graph_and_jobs = Load( $pipeline_graph_and_jobs_yaml );

print $pipeline_name, "\n";

#Get graph
#my $graph_easy = Graph::Easy->new( $graph_and_jobs->[1] );
#my $graph = Graph::Convert->as_graph ( $graph_easy );

my $graph = make_graph1($graph_and_jobs);
my @job_ids = $graph->vertices;
print "Got these jobs:  @job_ids\n";

#Get jobs
my @lsf_jobs;


#####################

# 1.
# Assign an lsf job name to each job - so we can refer to it in the dependency
# expression. This will be a combination of pipeline name, step, and a time based
# run number

my $timestamp = gmstamp();

foreach my $job_id (@job_ids) {
   $graph->set_vertex_attribute($job_id, 'lsf-job-name', $pipeline_name . "-$job_id-$timestamp" );
   print "Set $job_id lsf name to $pipeline_name-$job_id-$timestamp\n"
}

my $job_manager = LSF::JobManager->new( -q => 'normal', -M => "500" ,-R => "span[hosts=1] select[mem>500] rusage[mem=500]" , -L => "/usr/local/bin/bash");#, -L => '/bin/bash' );

# 2.
# Submit each job with it's possible dependencies.
my @job_ids_topologically_sorted = $graph->topological_sort;
foreach my $job_id ( @job_ids_topologically_sorted  ){
  print "Processing $job_id\n";
  my @dependent_jobs = $graph->predecessors($job_id);
  my $dependency_str;
  if( @dependent_jobs > 0 ) {
      print "Found dependent jobs for job $job_id: jobs(s) @dependent_jobs\n";
      my @dependencies;
      foreach my $dependent_job (@dependent_jobs ) {
         my $lsf_job_name = $graph->get_vertex_attribute( $dependent_job, 'lsf-job-name');
         push( @dependencies, "done($lsf_job_name)");
      }
      $dependency_str = join('&&', @dependencies);

      print "LSF dependency string: ", $dependency_str, "\n";
   }else{ print "No dependent job for $job_id.\n"}


   #Give the job a name with -J
   my $lsf_job_name = $graph->get_vertex_attribute( $job_id, 'lsf-job-name' );
   my ($job,$step) = split(/\./,$job_id);
   my $cmd = $graph_and_jobs->{ $job }->{$step}->{cmd};
   my $mem = $graph_and_jobs->{ $job }->{$step}->{mem};
   my $queue = $graph_and_jobs->{ $job }->{$step}->{queue};
   my $cores = $graph_and_jobs->{ $job }->{$step}->{cores};
   my %args;
   $args{-w} = $dependency_str if @dependent_jobs >0;
   $args{-J} = $lsf_job_name;
   $args{-R} = "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" if defined($mem);
   $args{-M} = $mem if defined($mem); # in KB, but we * by 1000 not 1024, since esub requires it to be less.
   $args{-o} = '/nfs/users/nfs_r/rg12/lustre110/Pipeline/lsfout/out-%J';
   $args{'-q'} = $queue if defined($queue);
   $args{'-n'} = $cores if defined($cores);
   print "\n---SUBMISSION---\n";
   if(  defined( $cmd ) ){
     print "CMD:\n$cmd\n";
   }else{
     print "CMD IS NOT DEFINED - (groupby or once step perhaps?)\n";
     $cmd = "echo nothing > stdout"
   }
   print "Submit $lsf_job_name with :\n" , join " ",  %args,"\n\n";
   $job_manager->submit(%args, $cmd);
   #if( @dependent_jobs > 0  ) {
   #  print "Submit $lsf_job_name\n";
   #  $job_manager->submit( -w => $dependency_str, -J => $lsf_job_name, $cmd);
   #} else {
   #  $job_manager->submit( -J => $lsf_job_name, $cmd);
   #}

}


#3. wait
for my $job ($job_manager->jobs){
    $job->top;
}

$job_manager->wait_all_children( history => 1 );
print "All children have completed!\n";

for my $job ($job_manager->jobs){ # much quicker if history is pre-cached
    print STDERR "$job exited non zero\n" if $job->history->exit_status != 0;

}

$job_manager->clear; # clear out the job manager to reuse.


sub make_graph1 {
   my $jobs = shift;
   my $graph = Graph->new;
   foreach my $job_num (%$jobs){
    my $stepsinjob = $jobs->{$job_num};
    foreach my $stepinjob ( keys %$stepsinjob ) {
       if( $stepsinjob->{$stepinjob}->{dependents}->length > 0 ){
          foreach my $dependency ($stepsinjob->{$stepinjob}->{dependents}->flatten) {
              $graph->add_edge( $dependency, "$job_num.$stepinjob" );
              print "added  $dependency - $job_num.$stepinjob  to graph\n";
          }
       }else{
        $graph->add_vertex("$job_num.$stepinjob");
        print "add vertex $job_num.$stepinjob\n";
       }
     }
   }
   return $graph;
}


sub make_graph {
   my $jobs = shift;
   my $graph = Graph->new;
   my $i = 0;
   foreach my $jobs (@$jobs){
     foreach my $stepinjob ( keys %$jobs ) {
       if( $jobs->{$stepinjob}->{dependents}->length > 0 ){
          foreach my $dependency ($jobs->{$stepinjob}->{dependents}->flatten) {
              $graph->add_edge( $dependency, "$i.$stepinjob" );
              print "added  $dependency - $i.$stepinjob  to graph\n";
          }
       }else{
        $graph->add_vertex("$i.$stepinjob");
        print "add vertex $i.$stepinjob\n";
       }
     }
     $i++;
   }
   return $graph;
}

