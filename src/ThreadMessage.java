/*
	This file is part of OSPREY.

	OSPREY Protein Redesign Software Version 2.1 beta
	Copyright (C) 2001-2012 Bruce Donald Lab, Duke University

	OSPREY is free software: you can redistribute it and/or modify
	it under the terms of the GNU Lesser General Public License as
	published by the Free Software Foundation, either version 3 of
	the License, or (at your option) any later version.

	OSPREY is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
	GNU Lesser General Public License for more details.

	You should have received a copy of the GNU Lesser General Public
	License along with this library; if not, see:
	      <http://www.gnu.org/licenses/>.

	There are additional restrictions imposed on the use and distribution
	of this open-source code, including: (A) this header must be included
	in any modification or extension of the code; (B) you are required to
	cite our papers in any publications that use this code. The citation
	for the various different modules of our software, together with a
	complete list of requirements and restrictions are found in the
	document license.pdf enclosed with this distribution.

	Contact Info:
			Bruce Donald
			Duke University
			Department of Computer Science
			Levine Science Research Center (LSRC)
			Durham
			NC 27708-0129
			USA
			e-mail:   www.cs.duke.edu/brd/

	<signature of Bruce Donald>, Mar 1, 2012
	Bruce Donald, Professor of Computer Science
*/
///////////////////////////////////////////////////////////////////////////////////////////////
//	ThreadMessage.java
//
//	Version:           2.1 beta
//
//
//	  authors:
// 	  initials    name                 organization                email
//	 ---------   -----------------    ------------------------    ----------------------------
//	  KER        Kyle E. Roberts       Duke University         ker17@duke.edu
///////////////////////////////////////////////////////////////////////////////////////////////
public class ThreadMessage {
    final static int UNDEFINED = -1 ;
    final static int NULL      =  0 ;
    final static int BYTE      =  1 ;
    final static int CHAR      =  2 ;
    final static int SHORT     =  3 ;
    final static int BOOLEAN   =  4 ;
    final static int INT       =  5 ;
    final static int LONG      =  6 ;
    final static int FLOAT     =  7 ;
    final static int DOUBLE    =  8 ;
    final static int PACKED    =  9 ;
    final static int LB        = 10 ;
    final static int UB        = 11 ;
    final static int OBJECT    = 12 ;

    final static int ANY_SOURCE = -1;
    final	static int ANY_TAG = -1;

    public int dataType;
    public int source;
    public int tag;
    public Object obj;




    public ThreadMessage(Object obj, int dataType, int source, int tag) {
        this.dataType = dataType;
        this.source = source;
        this.tag = tag;
        this.obj = obj;
    }




    public Object getObj(int offset, int count) {
        return determineObj(obj,dataType,offset,count);
    }

    public Object determineObj(Object obj, int dataType, int offset, int count) {
        Object o = null;
        switch(dataType) {
        case INT :
            o = new Integer[count];
            if(obj instanceof Integer[]) {
                for(int i=offset; i<offset+count; i++)
                    ((Integer[])o)[i-offset]= ((Integer[])obj)[i];
            } else {
                for(int i=offset; i<offset+count; i++)
                    ((Integer[])o)[i-offset]= ((int[])obj)[i];
            }
            return o;

        case BOOLEAN :
            o = new Boolean[count];
            if(obj instanceof Boolean[]) {
                for(int i=offset; i<offset+count; i++)
                    ((Boolean[])o)[i-offset]= ((Boolean[])obj)[i];
            } else {
                for(int i=offset; i<offset+count; i++)
                    ((Boolean[])o)[i-offset]= ((boolean[])obj)[i];
            }
            return o;
        case FLOAT :
            o = new Float[count];
            if(obj instanceof Float[]) {
                for(int i=offset; i<offset+count; i++)
                    ((Float[])o)[i-offset]= ((Float[])obj)[i];
            } else {
                for(int i=offset; i<offset+count; i++)
                    ((Float[])o)[i-offset]= ((float[])obj)[i];
            }
            return o;
        case DOUBLE :
            o = new double[count];
            for(int i=offset; i<offset+count; i++)
                ((double[])o)[i-offset]= ((double[])obj)[i];

            return o;
        case OBJECT :
            o = new Object[count];
            for(int i=offset; i<offset+count; i++)
                ((Object[])o)[i-offset]= ((Object[])obj)[i];

            return o;
        case BYTE :
            o = new byte[count];
            for(int i=offset; i<offset+count; i++)
                ((byte[])o)[i-offset]= ((byte[])obj)[i];

            return o;

        default:
            System.out.println("TRYING TO GET OBJ THAT HASN'T BEEN WRITTEN YET: "+dataType);
            System.exit(0);
        }

        return o;
    }


}
