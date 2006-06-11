<?
# -*- Mode: C++; tab-width: 2; -*-
# vi: set ts=2:
#
# --------------------------------------------------------------------------
#                   OpenMS Mass Spectrometry Framework
# --------------------------------------------------------------------------
#  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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
#
# --------------------------------------------------------------------------
# $Id: remove_include_guards.php,v 1.3 2005/04/06 09:36:57 marc_sturm Exp $
# $Author: marc_sturm $
# $Maintainer: Marc Sturm $
# --------------------------------------------------------------------------

include "common_functions.php";

	########################usage###############################
	if ($argc<2 OR $argc>3)
	{
?>
Usage: remove_include_guards.php <Path to OpenMS> [-v]

Removes all redundant include guards (those around include statements).

options:
  -v verbose output

<?
	exit;	
	}
	
	######################command line parsing###########################
	if (endsWith($argv[1],"/"))
	{
		$path = substr($argv[1],0,-1);
	}
	else
	{
		$path = $argv[1];
	}
	
	if ($argv[2] == "-v")
	{
		$verbose = true;
	}
	else
	{
		$verbose = false;
	}
	
	
	print "\n\n";
	#################include guards###############################
	$files=array();
	exec("find $path/source/ -name \"*.C\"", $files);
	exec("find $path/include/ -name \"*.h\"", $files); //todo
	foreach ($files as $header)
	{
		if ($verbose)
		{
			print "$header\n";
		}
		$out = array();
		$found = false;
		$file = file($header);
		for ($i=0;$i<count($file);$i++)
		{
			$line = trim($file[$i]);
			$nline = trim($file[$i+1]);
			$nnline = trim($file[$i+2]);
			if (beginsWith($line,"#ifndef OPENMS_") AND isIncludeLine($nline,$include) AND beginsWith($nnline,"#endif") )
			{
				$out[] = $file[$i+1];
				$i +=2;
				$found = true;
			}
			else
			{
				$out[] = $file[$i];
			}
		}
		
		if ($found)
		{
			print "WRITING: $header\n";
			writeFile($header,$out);
		}
	}
	



?>
