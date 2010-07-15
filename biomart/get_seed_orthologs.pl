#!/usr/bin/perl -w
#   Author: Francois Serra
# Creation Date: 2010/07/05 16:51:42





use strict;
use BioMart::Initializer;
use BioMart::Query;
use BioMart::QueryRunner;

my $confFile = "conf/metazoa_mart_reg.xml";
die ("Cant find configuration file $confFile\n") unless (-f $confFile);



my $action='cached';
my $initializer = BioMart::Initializer->new('registryFile'=>$confFile, 'action'=>$action);
my $registry = $initializer->getRegistry;

my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');

	$query->setDataset("dmelanogaster_eg_gene");
	$query->addFilter("status", ["KNOWN"]);
	$query->addAttribute("ensembl_gene_id");
	$query->addAttribute("ensembl_transcript_id");
	$query->addAttribute("coding");

$query->formatter("FASTA");

my $query_runner = BioMart::QueryRunner->new();
$query_runner->execute($query);
$query_runner->printHeader();
$query_runner->printResults();
$query_runner->printFooter();











#my $query = BioMart::Query->new('registry'=>$registry,'virtualSchemaName'=>'default');
#
#	$query->setDataset("dmelanogaster_eg_gene");
#	$query->addFilter("status", ["KNOWN"]);
#	$query->addAttribute("ensembl_gene_id");
#	$query->addAttribute("dananassae_eg_gene");
#	$query->addAttribute("derecta_eg_gene");
#
#$query->formatter("TSV");
#






























