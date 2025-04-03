#!/usr/bin/perl
foreach $ar (@ARGV)
{
    $r1 = $ar;
    $r2 = $ar;
    $r2 =~ s/_1/_2/;

    $name = $ar;
    $name =~ s/_1.fastq//;

    system("bowtie2 -p 20 -x /usr/local/BBDD/WoLr2/WoLr2 -1 $r1 -2 $r2 --very-sensitive --no-head --no-unal | cut -f1-9 | sed \"s/\$/\\t*\\t*/\" | gzip > $name.sam.gz");
}
