require 'fileutils'
require 'yaml'

ls_Version = '0.1.1'

puts 'Building GPF executables...'

ls_DestDir = "gpf-#{ls_Version}-linux"
FileUtils.rmtree(ls_DestDir)
FileUtils.rmtree('obj')
FileUtils.mkpath(ls_DestDir)

lk_Projects = ['gpfbatch', 'gpfd', 'gpfdump', 'gpfindex', 'gpfquery']
lk_Projects.each { |ls_Project| system("cd projects/#{ls_Project} && qmake && make clean && cd ../../") }
lk_Projects.each { |ls_Project| system("cd projects/#{ls_Project} && qmake && make release && cd ../../") }
lk_Projects.each { |ls_Project| FileUtils::cp(ls_Project, File::join(ls_DestDir, ls_Project)) }

system("tar cvf #{ls_DestDir}.tar #{ls_DestDir}")
system("bzip2 -z #{ls_DestDir}.tar")
