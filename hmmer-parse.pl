tie @go, 'Tie::File', $go_ref or die "Problem tying $go_ref";
#! /usr/bin/perl
use warnings;
use strict;
use Bio::SearchIO;
use Bio::Index::Fasta;

die "Usage: hmmer-parse.pl <hmmer output> <fasta file> <offset>/n" unless @ARGV==3;


my $file=shift @ARGV;
my $seqfile=shift @ARGV;
my $offset=shift @ARGV;

my $cutoff='20';
(my $name) = ($file =~ /([^\.]+)\.\w+/);



my $index = Bio::Index::Fasta->new('-filename' => $seqfile, '-write_flag' => 1);

my %Sequences=read_fasta($seqfile);

(my $sname) = ($seqfile =~ /([^\.]+)\.\w+/);


open(FAS,">$name\.hmm\.test\.fa");

my $count;

my $in = new Bio::SearchIO(-format => 'hmmer', 
						   -file   => $file);
while( my $result = $in->next_result ) {
        # this is a Bio::Search::Result::HMMERResult object
        print $result->query_name(), " for HMM ", $result->hmm_name(), "\n";
        while( my $hit = $result->next_hit ) {
            while( my $hsp = $hit->next_hsp ) {
            	if($hsp->score() > $cutoff) {
            		$count++; 
            		my($start,$end)=$hsp->range('query');
            		my $id=$hit->name();
            		my $score=$hsp->score();
            		my $seqstart=$start-$offset;
            		if ($seqstart < 0) {$seqstart=0; }
            		my $seq=$index->fetch($hit->name());
            		my $seqext=substr($seq->seq,$seqstart,($end-$start+2*$offset));
            		print FAS "\>$id\-$score\n$seqext\n";
				}
			}
        }
    }

print "$count hits !\n";

