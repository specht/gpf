require 'fileutils'
require 'yaml'

def determinePlatform()
	case RUBY_PLATFORM.downcase
	when /linux/
		'linux'
	when /darwin/
		'mac'
	when /mswin/
		'windows'
	else
		puts "Internal error: #{RUBY_PLATFORM} platform not supported."
		exit 1
	end
end


ls_Platform = determinePlatform()

ls_Config = ''

if ls_Platform == 'windows' 
	unless File::exists?('release.conf')
		puts 'Error: Need release.conf!'
		exit
	else
		ls_Config = File::read('release.conf')
	end
end

eval(ls_Config)

ls_Version = File::basename(Dir::pwd())

if ls_Version == nil
	puts 'Internal error: Unable to determine Proteomatic version.'
	exit 1
end

File.open('version.txt', 'w') { |lk_File| lk_File.write(ls_Version) }

puts "Creating release for GPF #{ls_Version}..."

ls_Make = {'windows' => 'make', 'linux' => 'make', 'mac' => 'make'}
ls_QMake = {'windows' => 'qmake -spec win32-g++', 'linux' => 'qmake', 'mac' => 'qmake -spec macx-g++'}
ls_BinaryExtension = {'windows' => '.exe', 'linux' => '', 'mac' => ''}

ls_DestDir = "gpf-#{ls_Version}-#{ls_Platform}"
FileUtils.rmtree(ls_DestDir)
FileUtils.mkpath(ls_DestDir)
FileUtils.rm_rf(ls_DestDir + '.tar')
FileUtils.rm_rf(ls_DestDir + '.tar.bz2')
FileUtils.rm_rf(ls_DestDir + '.zip')

puts 'Building GPF executables...'

FileUtils.rmtree(File::join('obj'))
lk_Projects = ['gpfbatch', 'gpfd', 'gpfdump', 'gpfindex', 'gpfquery']
lk_Projects.each { |ls_Project| system("cd projects/#{ls_Project} && #{ls_QMake[ls_Platform]} && #{ls_Make[ls_Platform]} release && cd ../../") }

puts 'Collecting GPF executables...'

lk_Projects.each { |ls_Project| FileUtils.cp(ls_Project + ls_BinaryExtension[ls_Platform], ls_DestDir) }

if (ls_Platform == 'windows')
	FileUtils.cp(File::join(QT_PATH, 'bin/QtCore4.dll'), ls_DestDir)
	FileUtils.cp(File::join(QT_PATH, 'bin/QtNetwork4.dll'), ls_DestDir)
	FileUtils.cp(File::join(MINGW_PATH, 'bin/mingwm10.dll'), ls_DestDir)
end

if (ls_Platform == 'windows')
	puts 'Building ZIP package...'
	
	system("helpers/7z/7za.exe a -r #{ls_DestDir}.zip #{ls_DestDir}")
else
	puts 'Builing bzip2 package...'
	system("tar cvf #{ls_DestDir}.tar #{ls_DestDir}")
	system("bzip2 -z #{ls_DestDir}.tar")
end

FileUtils.rmtree(ls_DestDir)
FileUtils.rmtree(File::join('obj'))
lk_Projects.each { |ls_Project| FileUtils::rm_rf(ls_Project + ls_BinaryExtension[ls_Platform]) }
