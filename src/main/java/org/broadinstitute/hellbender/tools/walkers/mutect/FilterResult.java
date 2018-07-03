package org.broadinstitute.hellbender.tools.walkers.mutect;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class FilterResult {
    private List<String> filtersApplied = new ArrayList<>();
    private Map<String, Object> newAnnotations = new HashMap<>();

    public void addFilters(String filter){
    }

}

