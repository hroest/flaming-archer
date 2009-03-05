<?php
 # -*- mode: php; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License, or (at your option) any later version.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
# --------------------------------------------------------------------------
# $Maintainer: Marc Sturm $
# $Authors: $
# --------------------------------------------------------------------------

///write array to file line by line
function writeFile($name,$array)
{
	$fp = fopen($name,"w");
	foreach ($array as $line){
		fwrite($fp,rtrim($line,"\r\n")."\n");
		}
	fclose($fp);
}

/// returns true if $string begins with $begin, return false otherwise
function beginsWith($string,$begin)
{
	if (substr($string,0,strlen($begin)) == $begin )
	{
		return true;
	}
	return false;
}

/// returns true if $string ends with $end, return false otherwise
function endsWith($string,$end)
{
	if (substr($string,-1*strlen($end)) == $end )
	{
		return true;
	}
	return false;
}

/// returns the prefix of length $length of string $string
function prefix($string,$length)
{
	return substr($string,$length);
}

/// returns the suffix of length $length of string $string
function suffix($string,$length)
{
	return substr($string,-1*$length);
}

/**
	returns the include guard for the included file $include
	give only the path and filename i.e. the part between "#include <" and ">"
*/
function includeToGuard($include)
{
	return strtoupper(str_replace(".","_",str_replace("/","_",$include)));
}

/**
	return true if the given line is a include line, returns false otherwise
	the fielname and path of the include are written into $include if it is given 
*/
function isIncludeLine($line,&$include)
{
	if (ereg ("^#[ \t]*include[ \t]*<(.*)>", ltrim($line),$parts))
	{
		$include = $parts[1];
		return true;
	}
	return false;
}

function isAlpha($char)
{
	$c = ord($char);
	return ((($c >= 64) && ($c <= 90)) || (($c >= 97) && ($c <= 122)));
}

function isSpace($char)
{
	return ($char == " ") || ($char == "\n") || ($char == "\t");
}

function isDigit($char)
{
	return ('0' <= $char) && ($char <= '9');
}

function isAlnum($char)
{
	return isAlpha($char) || isDigit($char);
}

function isIdentifierChar($char)
{
	return isAlnum($char) || $char == '_';
}

/// tokenize a line into something like C++-Tokens (not exactly, but good enough)
function tokenize($line)
{
	$result = array();
	$start = 0;
	// states: start => whitespace
	//         ident => identifier / number etc.
	$state = 0;
	for ($i = 0; $i != strlen($line); $i++)
	{
		$char = $line[$i];
		if ($state == 0)
		{
			// start state
			if (isIdentifierChar($char))
			{
				// identifier / number etc.
				$start = $i;
				$state = 1;
			}
			else if (!isSpace($char))
			{
				// operator
				$result[] = $char;
			}
		}
		else
		{
			// identifier state
			if (!isIdentifierChar($char))
			{
				// end of identifier
				$state = 0;
				$result[] = substr($line, $start, $i - $start);
				
				if (!isSpace($char))
				{
					// additionally, the current char belongs to an operator
					$result[] = $char;
				}
			}
		}
	}
	
	return $result;
}

/// parses out the maintainers from the maintainer line
function parseMaintainerLine($line)
{
	$replacements = array(
												 "Maintainer"=>"",
												 ":"=>"",
												 "/"=>"",
												 "$"=>"",
												 "	"=>" ",
												 "  "=>" ",
												 "#"=>""
												);
	$result = explode(",",strtr($line,$replacements));
	$result = array_map("trim",$result);
	$result = array_filter($result,"strlen");
	return array_unique($result);
}

/**
	@brief parses the member information form a file
	
	@param path The OpenMS path
	@param header The header file name
	
	@return Returns the parses information
*/
function getClassInfo($bin_path,$header, $debug)
{
	$members = array(
		"classname" => substr(basename($header),0,-2),
		"public-long" => array(),
		"public" => array(),
		"non-public" => array(),
		"variables" => array(),
		);

	######################## needed stuff ###############################
	if (!file_exists("$bin_path/doc/xml_output/"))
	{
		print "Error: The directory '$bin_path/doc/xml_output/' is needed!\n";
		print "       Please execute 'make idoc' first!\n";
	}
	
	######################## load file ###############################
	$paths = array(
		"",
		"Internal_1_1",
		"Math_1_1",
		);
	
	$found = false;
	foreach ($paths as $p)
	{
		//find class
		$tmp = "$bin_path/doc/xml_output/classOpenMS_1_1".$p.$members["classname"].".xml";
		if (file_exists($tmp))
		{
			$class = simplexml_load_file($tmp);
			$found = true;
			break;
		}
		//find struct
		else
		{
			$tmp = "$bin_path/doc/xml_output/structOpenMS_1_1".$p.$members["classname"].".xml";
			if (file_exists($tmp))
			{
				$class = simplexml_load_file($tmp);
				$found = true;
				break;
			}
			
		}
	}
	if (!$found)
	{
		print "Error: No XML file found for class '".$members["classname"]."'. Aborting!\n";
		return $members;
	}
	
	######################## parse ###############################
	
	$classname = substr($class->compounddef->compoundname,8);
	foreach ($class->compounddef->sectiondef as $section)
	{
		foreach($section->memberdef as $member)
		{
			#method
			if ($member["kind"]=="function")
			{
				#public methods
				if ($member["prot"]=="public")
				{
					#name
					$mem = $member->definition.$member->argsstring;
					
					#template parameters
					if (isset($member->templateparamlist))
					{
						$first = true;
						$template = "";
						foreach($member->templateparamlist->param as $para)
						{
							# edited by Clemens 2009-01-15; this seems to work for doxygen 1.5.7.1
							# (handling of whitespace is horrific, please improve if you know better)
							if ($first)
							{
								$template = "template <";
								$first = false;
							}
							else
							{
								$template .= ",";
							}
							if (isset($para->type))
							{
								$template .= " ".$para->type;
								$template = trim($template);
							}
							if (isset($para->type->ref))
							{
								$template .= " ".$para->type->ref;
								$template = trim($template);
							}
							if (isset($para->defname))
							{
								$template .= " ".$para->defname;
								$template = trim($template);
							}
						}
						if ( !$first )
						{
							$template .= " > ";
						}
						$mem = $template.$mem;
					}
					
					if ($debug>4)
					{
						print "Exceptions for '$mem':\n";
						var_dump($member->exceptions);
					}
					
					#exceptions
					if (isset($member->exceptions))
					{
						$except = " throw (";
						$first = true;
						foreach($member->exceptions->ref as $ref)
						{
							if (!$first)
							{
								$except .= ", ";
							}
							else
							{
								$first = false;
							}
							$except .= (string) $ref;
						}
						if ($except!=" throw (")
						{
							$mem .= $except.")";
						}
					}
					# remove namespace stuff
					$mem = strtr($mem,array($classname."::"=>""));
					
					$members["public-long"][]=$mem;
					$members["public"][] = (string)$member->name;
				}
				#non-public methods
				else
				{
					$members["non-public"][] = (string)$member->name;
				}
			}
			else if($member["kind"]=="variable" && $member["prot"]!="public")
			{
				$members["variables"][] = (string)$member->name;
			}
		}
	}
	return $members;
}

/// Load information about the tested methods
function parseTestFile($filename)
{
	$tests = array();
	$todo_tests = array();
	$in_test = false;
	$last_test = "";
	#load file
	$tmp = file($filename);
	
	foreach ($tmp as $line)
	{	
		$line = trim($line);
		if (beginsWith($line,"START_SECTION(") || beginsWith($line,"START_SECTION (") || beginsWith($line,"START_SECTION	("))
		{
			# strip brackets
			$function = trim(substr($line,13));
			//print "Function: $function \n";
			while ($function[0]=='(' && $function[strlen($function)-1]==')')
			{
				$function = trim(substr($function,1,-1));
			}
			# ignore extra function tests
			if (!beginsWith($function,"[EXTRA]"))
			{
				$tests[] = $function;
				$last_test = $function;
			}
			$in_test = true;
		}
		else if (beginsWith($line,"END_SECTION"))
		{
			$in_test = false;
		}
		else if ($in_test)
		{
			if (strpos($line,"TODO")!==FALSE || strpos($line,"????")!==FALSE)
			{
				$todo_tests[] = $function;
			}
		}
	}
	return array("todo"=>$todo_tests,"tests"=>$tests);
}

/// Comares declared and tested methods
function compareDeclarationsAndTests($declarations,$tests)
{
	# Replacements to use when comparing methods
	$method_replacements = array(
	"std::" => "",
	"OpenMS::" => "",
	" " => "",
	"\t" => "",
	"=0" => "",
	"throw()" => "",
	"virtual" => "",
	"static" => "",
	"/*" => "",
	"*/" => "",
	);
	
	$done = array();
	
	$out = array(
		"missing" => array(),
		"unknown" => array(),
		"double"  => array()
		);

	#make a copy without whitespaces
	$tmp  = array();
	foreach ($declarations as $m)
	{
		$tmp[] = strtr($m,$method_replacements);
	}
	
	#compare tests and declarations
	foreach ($tests as $t)
	{
		$stripped = strtr($t,$method_replacements);
		$pos = array_search($stripped,$tmp);
		if ($pos === FALSE)
		{
			if (in_array($stripped,$done))
			{
				$out["double"][] = $t;
			}
			else
			{
				$out["unknown"][] = $t;
			}
		}
		else
		{
			unset($tmp[$pos]);
			$done[] = $stripped;
		}
	}
	
	#report extra tests
	if (count($tmp)!=0)
	{
		foreach ($tmp as $t)
		{
			# look up test with spaces
			foreach($declarations as $z)
			{
				if (strtr($z,$method_replacements)==$t)
				{
					$out["missing"][] = $z;	
					break;
				}
			}
		}
	}
	return $out;
}

?>
