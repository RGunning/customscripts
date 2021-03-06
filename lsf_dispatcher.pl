#!/usr/bin/env perl
#===============================================================================
#
#					FILE: lsf_job_manager.pl
#
#				 USAGE: ./lsf_job_manager.pl
#
#	 DESCRIPTION:
#
#			 OPTIONS: ---
# REQUIREMENTS: ---
#					BUGS: ---
#				 NOTES: ---
#				AUTHOR: NJWALKER (), nw11@sanger.ac.uk modified RJGUNNING, rg12@sanger.ac.uk
#			 COMPANY:
#			 VERSION: 1.0
#			 CREATED: 24/02/13 12:05:44
#			REVISION: ---
#===============================================================================

use Modern::Perl;
use YAML::Any;
use Graph;
use Path::Tiny;
use Data::Dumper;


#use Graph::Easy;
#use Graph;
#use Graph::Convert;
#use Path::Class qw(dir file);
#use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
#use LSF::Job;
#use autobox::Core;

# BASIC BLOCKING DISPATCHER
#
# Currently read yaml file into graph - then go through the job_ids and see
# if they have a dependency - if not execute and delete from graph - if do skip
# and keep going until no graph left. Probably better way to do this breadth first
# and focus on leaves

die "Did not receive plite yaml file " unless defined($ARGV[0]);
my $plite_yaml = path($ARGV[0]);
my $pipeline_graph_and_jobs_yaml = $plite_yaml->slurp;
my $graph_and_jobs = Load( $pipeline_graph_and_jobs_yaml );
my $pipeline_name = $plite_yaml->basename('.graph.yaml');

my $graph = make_graph( $graph_and_jobs);
#my @job_ids = $graph->vertices;
my @job_ids = $graph->topological_sort;
dispatch_lsf(\@job_ids, $graph, $pipeline_name);

sub make_graph {
	my $jobs = shift;
	my $graph = Graph->new;

	foreach my $job_num (keys %$jobs){
		my $job = $jobs->{$job_num};
		foreach my $stepinjob ( keys %$job ) {
			my $dependents = $job->{$stepinjob}->{dependents};
			if( @$dependents > 0 ){
				foreach my $dependency (@$dependents) {
					$graph->add_edge( $dependency, "$job_num.$stepinjob" );	 #print "added	$dependency - $job_num.$stepinjob	 to graph\n";
				}
			}else{
				$graph->add_vertex("$job_num.$stepinjob"); #print "add vertex $job_num.$stepinjob\n";
			}
		}
	}
	return $graph;
}

sub dispatch_lsf {
	my ($job_ids, $graph, $pipeline_name ) = @_;
	my $job_id;

	use Time::Stamp qw(-stamps);
	use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;

	my $job_manager = LSF::JobManager->new( -q => 'normal', -M => "500" , -n => "1" , -R => "span[hosts=1] select[mem>500] rusage[mem=500]" , -L => "/usr/local/bin/bash");#, -L => '/bin/bash' );

	my $timestamp = gmstamp();
	foreach my $job_id (@$job_ids) {
		$graph->set_vertex_attribute($job_id, 'lsf-job-name', $pipeline_name . "-$job_id-$timestamp" );
		print STDERR "Set $job_id lsf name to $pipeline_name-$job_id-$timestamp\n"
	}

	while( @$job_ids > 0 ){
		$job_id = shift @$job_ids;

		my @dependent_jobs = $graph->predecessors($job_id);
		my $dependency_str;

		if( @dependent_jobs > 0 ) {
			print STDERR "Found dependent jobs for job $job_id: jobs(s) @dependent_jobs\n";
			my @dependencies;
			foreach my $dependent_job (@dependent_jobs ) {
				my $lsf_job_name = $graph->get_vertex_attribute( $dependent_job, 'lsf-job-name');
				push( @dependencies, "done($lsf_job_name)");
			}
			$dependency_str = join('&&', @dependencies);

			print "LSF dependency string: ", $dependency_str, "\n";
	 	} else { print "No dependent job for $job_id.\n"}

		#Give the job a name with -J
		my $lsf_job_name = $graph->get_vertex_attribute( $job_id, 'lsf-job-name' );
		my ($job,$step) = split(/\./,$job_id);
	 	my ($cmd,$errorfile) = split(/2>/,$graph_and_jobs->{ $job }->{$step}->{cmd});
   		my $mem = $graph_and_jobs->{ $job }->{$step}->{mem};
   		my $queue = $graph_and_jobs->{ $job }->{$step}->{queue};
   		my $cores = $graph_and_jobs->{ $job }->{$step}->{cores};
   		my %args;
   		$args{-e} = $errorfile if defined($errorfile);
	 	$args{-w} = $dependency_str if @dependent_jobs >0;
	 	$args{-J} = $lsf_job_name;
		$args{-R} = "span[hosts=1] select[mem>$mem] rusage[mem=$mem]" if defined($mem);
		$args{-M} = $mem if defined($mem); # in KB, but we * by 1000 not 1024, since esub requires it to be less.
		$args{-o} = '/nfs/users/nfs_r/rg12/lustre110/Pipeline/lsfout/out-%J';
		$args{'-q'} = $queue if defined($queue);
		$args{'-n'} = $cores if defined($cores);
		print "\n---SUBMISSION---\n";
		if(	defined( $cmd ) ){
			print "CMD:\n$cmd\n";
		}else{
			print "CMD IS NOT DEFINED - (groupby or once step perhaps?)\n";
			$cmd = "echo nothing > stdout"
		}
		print "Submit $lsf_job_name with :\n" , join " ",	%args,"\n\n";
		$job_manager->submit(%args, $cmd);
	}
}

