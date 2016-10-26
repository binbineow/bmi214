

#!/usr/bin/env python

class SNP():
    '''class for chrY SNP queries in 1000G data.
       instructions:
         1. construct a SNP object with site info.
         2. walk vcf file:
            a. set ref/alt alleles for this site with setAlleles()
            b. track ref/alt/missing alleles across samples
    ''' 
    def __init__(self, snpInfoList):
        self.info = snpInfoList
        self.numMiss, self.numRef, self.numAlt = 0, 0, 0
        self.sampleList = []
        self.ref = None
    def setAlleles(self, ref, alt):
        self.ref = ref
        self.alt = alt
    def seeMiss(self): 
        'increment missingness count'
        self.numMiss += 1
    def seeRef(self): 
        'increment ref allele count'
        self.numRef += 1
    def seeAlt(self, sampleInfoTuple): 
        'record observed alt allele'
        self.numAlt += 1
        self.sampleList.append(sampleInfoTuple)
    def __str__(self):
        report = 'SNP: ' + ' '.join(self.info) + '\n'
        if self.ref is None:
            report += 'not found\n'
        else:
            report += '%3d missing\n'        % self.numMiss + \
                      '%3d %s reference\n'   % (self.numRef, self.ref) + \
                      '%3d %s alternative\n' % (self.numAlt, self.alt)

            report += 'Observed Derived Alleles In:\n'
            for sample in self.sampleList:
                report += '%s %s %s\n' % sample
        return report 

def main():
    'a contrived example to show how the SNP class works'

    site = SNP(['M9', 'K', 'rs3900', '20189645', '21730257'])
    site.setAlleles('C', 'G') # normally, this would be a dictionary look-up
    for _ in xrange(8): # pretend we observe 8 reference alleles
        site.seeRef()
    for _ in xrange(2): # pretend we observe 2 missing alleles
        site.seeMiss()          
    site.seeAlt(('NA19076', 'JPT', 'O3a2')) # observe a missing allele
    site.seeAlt(('NA19655', 'MXL', 'T1a1b'))
    print site

if __name__ == '__main__':
    main()