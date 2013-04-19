def db_schema
  ActiveRecord::Schema.define do
    unless table_exists? :strains
      create_table :strains do |t|
        t.column :name, :string
        t.column :description, :string
      end
    end

    unless table_exists? :features_snps
      create_table :features_snps, :id => false do |t|
        t.column :feature_id, :integer, :null => false
        t.column :snp_id, :integer, :null => false
      end
    end

    unless table_exists? :features
      create_table :features do |t|
        t.column :name, :string
        t.column :start, :integer
        t.column :end, :integer
        t.column :strand, :integer
        t.column :sequence, :string
      end
    end

    unless table_exists? :snps
      create_table :snps do |t|
        t.column :ref_pos, :integer
        t.column :qual, :float 
        t.column :reference_allele_id, :integer
      end
    end

    unless table_exists? :alleles
      create_table :alleles do |t|
        t.column :snp_id, :integer
        t.column :base, :string
      end
    end
   
    unless table_exists? :genotypes
      create_table :genotypes do |t|
        t.column :allele_id, :integer
        t.column :strain_id, :integer
        t.column :geno_qual, :float
      end
    end

    unless table_exists? :annotations
      create_table :annotations do |t|
        t.column  :qualifier, :string
        t.column :value, :string
        t.column :feature_id, :integer
      end
    end
     
    # indices
    unless index_exists? :strains, :id
      add_index :strains, :id
    end
    unless index_exists? :strains, :name
      add_index :strains, :name
    end
    unless index_exists? :features, :id
      add_index :features, :id
    end
    unless index_exists? :features, :name
      add_index :features, :name
    end
    unless index_exists? :features, :start
      add_index :features, :start
    end
    unless index_exists? :features, :end
      add_index :features, :end
    end
    unless index_exists? :features, :strand
      add_index :features, :strand
    end
    unless index_exists? :features, :sequence
      add_index :features, :sequence
    end
    unless index_exists? :snps, :id
      add_index :snps, :id
    end
    unless index_exists? :snps, :ref_pos
      add_index :snps, :ref_pos
    end
    unless index_exists? :snps, :qual
      add_index :snps, :qual
    end
    unless index_exists? :snps, :reference_allele_id
      add_index :snps, :reference_allele_id
    end
    unless index_exists? :features_snps, :feature_id
      add_index :features_snps, :feature_id
    end
    unless index_exists? :features_snps, :snp_id
      add_index :features_snps, :snp_id
    end
    unless index_exists? :snps, :qual
      add_index :snps, :qual
    end
    unless index_exists? :alleles, :snp_id
      add_index :alleles, :snp_id
    end
    unless index_exists? :alleles, :base
      add_index :alleles, :base
    end
    unless index_exists? :genotypes, :id
      add_index :genotypes, :id
    end
    unless index_exists? :genotypes, :allele_id
      add_index :genotypes, :allele_id
    end
    unless index_exists? :genotypes, :strain_id
      add_index :genotypes, :strain_id
    end
    unless index_exists? :genotypes, :geno_qual
      add_index :genotypes, :geno_qual
    end
    unless index_exists? :annotations, :feature_id
        add_index :annotations, :feature_id
    end
    unless index_exists? :annotations, :qualifier
      add_index :annotations, :qualifier
    end
    unless index_exists? :annotations, :value
      add_index :annotations, :value
    end
  end
end