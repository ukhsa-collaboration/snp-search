require 'snp_db_connection'

class Strain < ActiveRecord::Base
  has_many :alleles, :through => :genotypes
  has_many :genotypes
end

class Feature < ActiveRecord::Base
  has_and_belongs_to_many :snps
  has_many :annotations
end

class Snp < ActiveRecord::Base
  has_and_belongs_to_many :features
  has_many :alleles
  belongs_to :reference_allele, :class_name => "Allele", :foreign_key => "reference_allele_id"
end

class Allele < ActiveRecord::Base
  has_many :genotypes
  belongs_to :snp
  has_many :strains, :through => :genotypes
end

class Genotype < ActiveRecord::Base
  belongs_to :allele
  belongs_to :strain
end

class Annotation < ActiveRecord::Base
  belongs_to :feature
end
