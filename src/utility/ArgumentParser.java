package utility;

// from http://www.devx.com/tips/Tip/13004

//TODO: i don't like having to associate key-value pairs with =

import java.util.*;
public class ArgumentParser {
	
    private Map<String,List<String>> options = new Hashtable<String,List<String>>();
    
    public ArgumentParser(String[] args) {
    	List<String> arguments = null;
    	String key = null;
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("--")) {
            	if (key != null)
            	options.put(key.toLowerCase(), arguments);
                	key = args[i].substring(2);
                arguments = new LinkedList<String>();
            } else {
            	if (arguments == null) {
            		System.out.println("ERROR: malformed command-line input");
            		System.exit(-1);
            	}
                arguments.add(args[i]);
            }
        }
        if (key != null)
        	options.put(key.toLowerCase(), arguments);
    }

    public boolean hasOption(String opt) {
        return options.containsKey(opt.toLowerCase());
    }
    
    public boolean hasArguments(String opt) {
    	if (!hasOption(opt))
    		return false;
    	return (options.get(opt.toLowerCase()).size() > 0);
    }

    public List<String> getArguments(String opt) {
        return options.get(opt.toLowerCase());
    }
    
    public List<Pair<String, String[]>> getCLParserData() {
    	List<Pair<String, String[]>> result = new ArrayList<Pair<String,String[]>>();
    	
    	for (String s : options.keySet()) {
			int numArgs = options.get(s).size();
    		String args[] = new String[numArgs];
    		for (int i = 0; i < numArgs; i++)
    			args[i] = options.get(s).get(i);
    		result.add(new Pair<String,String[]>("--" + s, args));
    	}
    	
    	return result;
    }
    
    public String getArgument(String opt) {
    	if (hasArguments(opt))
    		return options.get(opt.toLowerCase()).get(0);
    	else
    		return null;
    }

}
