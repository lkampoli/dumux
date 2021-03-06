#!/bin/sh

# check if help is needed
if test "$1" = "--help" || test "$1" = "-help" \
   || test "$1" = "help" || test "$1" = ""; then
  echo ""
  echo "USAGE: $0 FILENAME"
  echo ""
  echo "The argument should be a C++ header or source file."
  echo "The file is converted to Markdown and forwarded to stdout."
  echo "In particular, a C++ comment is converted to its content,"
  echo "while normal code is put into corresponding Markdown code blocks."
  echo "Supposed to be called by the script merge_cpp_and_md.sh."
  exit 0
fi

if [[ ${1: -3} == ".hh" ]]; then
  sed '1,/#define/d' $1 >tmpfile
  isheader=true
else
  sed '1,/\*\*\*\*\*\//d' $1 >tmpfile
  isheader=false
fi

insidecodeblock=false
lastline=""
strippedlastline=""
firstline=true

while IFS= read -r line
do
  strippedline=$(echo $line | sed "s/^[ \t]*//")
  if [[ $firstline == false ]]; then
    if [[ ${strippedline:0:2} == "//" ]]; then
      if [[ $insidecodeblock == true ]]; then
        if [[ $lastline != "" ]]; then
          echo "$lastline"
        fi
        echo "\`\`\`"
        insidecodeblock=false
      else
        linetoprint=$(echo ${strippedlastline:2} | sed "s/^[ \t]*//")
        echo "$linetoprint"
      fi
    else
      if [[ $insidecodeblock == false  ]]; then
        linetoprint=$(echo ${strippedlastline:2} | sed "s/^[ \t]*//")
        echo "$linetoprint"
        echo "\`\`\`cpp"
        insidecodeblock=true
      else
        echo "$lastline"
      fi
    fi
  else
    if [[ ${strippedline:0:2} != "//" && $line != "" ]]; then
      echo "\`\`\`cpp"
      insidecodeblock=true
    fi
    firstline=false
  fi
  lastline="$line"
  strippedlastline="$strippedline"
done < tmpfile

if [[ $isheader == false ]]; then
  if [[ ${strippedlastline:0:2} == "//" ]]; then
    if [[ $insidecodeblock == true ]]; then
      echo "\`\`\`"
      insidecodeblock=false
    fi
    linetoprint=$(echo ${strippedlastline:2} | sed "s/^[ \t]*//")
    echo "$linetoprint"
  else
    if [[ $insidecodeblock == false  ]]; then
      echo "\`\`\`cpp"
      insidecodeblock=true
    fi
    echo "$lastline"
  fi
fi

if [[ $insidecodeblock == true ]]; then
  echo "\`\`\`"
  insidecodeblock=false
fi
