/*
 * File         : User_function.java
 * Author       : Xinying Li
 * Created      : Thu Apr  9 12:22:15 1998
 * Revision     : $Revision: 1.1 $
 * Updated      : $Date: 2007/12/11 18:44:48 $
 * Copyright: Northeast Parallel Architectures Center
 *            at Syracuse University 1998
 */

package mpi;

public abstract class User_function {

    /**
     * User-defined function for a new <tt>Op</tt>.
     * <p>
     * <table>
     * <tr><td><tt> invec       </tt></td><td> array of values to combine with
     *                                         <tt>inoutvec</tt> elements </tr>
     * <tr><td><tt> inoffset    </tt></td><td> initial offset in
     *                                         <tt>invec<tt> </tr>
     * <tr><td><tt> inoutvec    </tt></td><td> in-out array of accumulator
     *                                         locations </tr>
     * <tr><td><tt> inoutoffset </tt></td><td> initial offset in
     *                                         <tt>inoutvec<tt> </tr>
     * <tr><td><tt> count       </tt></td><td> number of items in arrays </tr>
     * <tr><td><tt> datatype    </tt></td><td> type of each item </tr>
     * </table>
     * <p>
     * Java equivalent of the MPI <tt>USER_FUNCTION</tt>.
     */

    public abstract void Call(Object invec, int inoffset,
                              Object inoutvec, int inoutoffset,
                              int count, Datatype datatype) ;
}

