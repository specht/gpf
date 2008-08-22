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
		puts 'Using default config. Please adjust values in release.conf if neccessary.'
		File::open('release.conf', 'w') do |lk_File|
			lk_File.puts "QT_PATH = 'c:/Qt/4.4.1'"
			lk_File.puts "MINGW_PATH = 'c:/MinGW'"
		end
	end
	ls_Config = File::read('release.conf')
end

eval(ls_Config)

# test config
if ls_Platform == 'windows' 
	lk_Errors = Array.new
	lk_Errors.push("Unable to find MinGW in #{MINGW_PATH}.") unless File::exists?(File::join(MINGW_PATH, 'bin/mingwm10.dll'))
	lk_Errors.push("Unable to find Qt in #{QT_PATH}.") unless File::exists?(File::join(QT_PATH, 'bin/QtCore4.dll'))
	unless lk_Errors.empty?
		puts 'Errors:'
		puts lk_Errors.join("\n")
		exit 1
	end
end

ls_Version = File::basename(Dir::pwd())

if ls_Version == nil
	puts 'Internal error: Unable to determine Proteomatic version.'
	exit 1
end

File.open('version.txt', 'w') { |lk_File| lk_File.write(ls_Version) }

puts "Creating release for GPF #{ls_Version}..."

ls_Make = {'windows' => File::join(MINGW_PATH, 'make.exe'), 'linux' => 'make', 'mac' => 'make'}
ls_QMake = {'windows' => File::join(QT_PATH, 'qmake') + ' -spec win32-g++', 'linux' => 'qmake', 'mac' => 'qmake -spec macx-g++'}
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
