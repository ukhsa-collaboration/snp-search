gem 'activerecord', "~> 3.1.3" 
require 'active_record'
gem 'sqlite3', "~> 1.3.4"
require 'sqlite3'
def establish_connection(db_location)
  ActiveRecord::Base.establish_connection(
    :adapter => "sqlite3",
    :database => db_location)
end