require 'rubygems'
require 'bio'
require  'snp_db_models'
require 'activerecord-import'
require 'diff/lcs'
require 'create_methods'
require 'filter_ignore_snps_methods'
require 'output_information_methods'

def find_unqiue_snps(strain_names, out, cuttoff_genotype, cuttoff_snp)

  *strain_names = strain_names

  where_statement = strain_names.collect{|strain_name| "strains.name = '#{strain_name}' OR "}.join("").sub(/ OR $/, "")

  outfile = File.open(out, "w")

   snps = Snp.find_by_sql("SELECT snps.* from snps INNER JOIN alleles ON alleles.snp_id = snps.id INNER JOIN genotypes ON alleles.id = genotypes.allele_id INNER JOIN strains ON strains.id = genotypes.strain_id WHERE (#{where_statement}) AND alleles.id <> snps.reference_allele_id AND genotypes.geno_qual >= #{cuttoff_genotype} AND snps.qual >= #{cuttoff_snp} AND (SELECT COUNT(*) from snps AS s INNER JOIN alleles ON alleles.snp_id = snps.id INNER JOIN genotypes ON alleles.id = genotypes.allele_id WHERE alleles.id <> snps.reference_allele_id and s.id = snps.id) = #{strain_names.size} GROUP BY snps.id HAVING COUNT(*) = #{strain_names.size}")
   # puts "The number of unique snps are #{snps.size}"

   output_information_methods(snps, outfile, cuttoff_genotype, cuttoff_snp, false)
end


def information(out, cuttoff_genotype, cuttoff_snp)

  puts "outputting SNP info....."
  
  strains = Strain.all

  snps = Snp.find_by_sql("SELECT distinct snps.* from snps INNER JOIN alleles ON alleles.snp_id = snps.id INNER JOIN genotypes ON alleles.id = genotypes.allele_id INNER JOIN strains ON strains.id = genotypes.strain_id where alleles.id <> snps.reference_allele_id")

  outfile = File.open(out, "w")

  output_information_methods(snps, outfile, cuttoff_genotype, cuttoff_snp, true)

end

