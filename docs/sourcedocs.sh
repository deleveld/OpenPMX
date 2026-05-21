#!/bin/sh
########################################################################
#
# This file is part of OpenPMX (https://github.com/deleveld/openpmx).
# Copyright (c) 2024 Douglas Eleveld.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

set -e

extract_md() 
{
	dir=${1}
	head=${2}
	
	echo '# ' $head
	echo
	find $dir -type f | sort |
	while IFS= read -r f; do
		awk -v file="$f" '
			BEGIN {
				sub(/^(\.\.\/|\.\/)+/, "", file)
			}
			function is_list_item(line) {
				return line ~ /^[ \t]*([*+-]|[0-9]+\.)[ \t]+/
			}
			function flush() {
				if (buf) {
					out = out buf
					if (!in_list)
						out = out "\n"
					buf = ""
				}
			}
			/^[[:space:]]*\/\/\// {
				line = $0
				sub(/^[[:space:]]*\/\/\/[ ]?/, "", line)

				if (is_list_item(line))
					in_list = 1
				else if (line !~ /^[ \t]/ && line != "")
					in_list = 0

				buf = buf line "\n"
				found = 1
				next
			}
			{ flush() }
			END {
				flush()
				if (found)
					printf "\n# %s\n\n%s", file, out
			}' "$f"
	done
}

extract_md ../include "" >SOURCEDOCS.md
extract_md ../src "" >>SOURCEDOCS.md
