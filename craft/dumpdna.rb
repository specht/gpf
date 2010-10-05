# Copyright (c) 2007-2010 Michael Specht
# 
# This file is part of GPF.
# 
# GPF is free software: you can redistribute it and/or modify it under 
# the terms of the GNU Lesser General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your 
# option) any later version.
# 
# GPF is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public 
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public 
# License along with GPF.  If not, see <http://www.gnu.org/licenses/>.

require 'yaml'

tripletToAminoAcid = Hash.new

File::open('dna-to-amino-acid.csv', 'r') do |f|
    f.each_line do |line|
        lineArray = line.split(';')
        triplet = lineArray.first.strip
        aa = lineArray.last.strip
        tripletToAminoAcid[triplet] = aa
    end
end

entries = Array.new
id = nil
entry = ''

File::open(ARGV.first, 'r') do |f|
    f.each_line do |line|
        line.strip!
        if line[0, 1] == '>'
            entries << [id, entry] if id && (!entry.empty?)
            id = line[1, line.size - 1]
            entry = ''
        else
            entry += line
        end
    end
end
entries << [id, entry] if id && (!entry.empty?)

width = 120

File::open('out.html', 'w') do |f|
    f.puts "<html>"
    f.puts "<head>"
    f.puts "</head>"
    f.puts "<body style='font-size: 8pt;'>"
    entries.each do |x|
        id = x[0]
        sequence = x[1]
        f.puts "<h1>#{id}</h1>"
        f.puts "<pre>"
        frames = [' ', ' ', ' ', ' ', ' ', ' ']
        i = 0
        currentFrame = 0
        while i <= sequence.size - 3
            triplet = sequence[i, 3].gsub('T', 'U')
            aa = tripletToAminoAcid[triplet]
            frames[currentFrame] += aa
            frames[(currentFrame + 1) % 3] += ' '
            frames[(currentFrame + 2) % 3] += ' '
            
            triplet = sequence[i, 3].gsub('T', 'U').reverse.tr('ACGU', 'UGCA')
            aa = tripletToAminoAcid[triplet]
            frames[((5 - currentFrame) % 3) + 3] += aa
            frames[((((5 - currentFrame) % 3) + 1) % 3) + 3] += ' '
            frames[((((5 - currentFrame) % 3) + 2) % 3) + 3] += ' '

            i += 1
            currentFrame = (currentFrame + 1) % 3
        end
        offset = 0
        while (!sequence.empty?)
            width.times { f.print('-') }
            f.puts
            thisWidth = 0
            while (thisWidth < width)
                num = "#{offset + thisWidth}"
                f.print(num)
                (5 - num.size).times { f.print(' ') }
                f.print('.    ')
                thisWidth += 10
            end
            f.puts
            sequencePart = sequence.slice!(0, width)
            (0...sequencePart.size).each do |index|
                f.print(sequencePart[index, 1])
            end
            f.puts
            (0...sequencePart.size).each do |index|
                f.print(sequencePart[index, 1].tr('ACGT', 'TGCA'))
            end
            f.puts
            (0..5).each { |x| f.puts frames[x].slice!(0, width) }
            offset += width
        end
        f.puts "</pre>"
    end
    f.puts "</body>"
    f.puts "</html>"
end
