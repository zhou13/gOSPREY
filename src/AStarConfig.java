public class AStarConfig {
    public static boolean enableJava	= true;
    public static boolean enableNativeC	= false;
    public static boolean enableCUDA	= false;
    public static long maxCPUMemory		= 1024L * 1024L * 1024L;  // 1G
    public static long maxGPUMemory		= 1024L * 1024L * 1024L;  // 1G

    public static int numWorkGroup		= 4;
    public static int numWorkItem		= 192;
    public static int numWorkItem2		= 192;

    public static double shrinkRatio	= 1.f;
};
