# plato-correlate

> Cross-correlation and dependency detection between PLATO tile streams

## What This Does

plato-correlate detects relationships between tile streams using cross-correlation. Given multiple named time series, it finds which ones are correlated, at what lag, and builds a dependency graph showing which streams predict which others.

## The Key Idea

If kitchen temperature rises 5 minutes before the HVAC kicks on, that's a lagged correlation. plato-correlate slides one time series past another at different lags and computes the Pearson correlation at each position. The lag with the strongest correlation tells you both the relationship strength and its timing. From all pairwise correlations, it builds a dependency graph: directed edges showing which streams predict which.

## Install

```bash
cargo add plato-correlate
```

## Quick Start

```rust
use plato_correlate::*;
use std::collections::HashMap;

let mut series = HashMap::new();
series.insert("temp".into(), vec![20.0, 21.0, 22.0, 23.0, 24.0]);
series.insert("hvac".into(), vec![0.0, 0.0, 1.0, 1.0, 1.0]);

// Build dependency graph with threshold
let graph = DependencyGraph::from_correlations(&series, 0.5);
for edge in &graph.edges {
    println!("{} → {} (strength={:.2}, lag={})", 
        edge.from, edge.to, edge.strength, edge.lag);
}

// Cross-correlate two specific series
let cc = cross_correlation(&series["temp"], &series["hvac"], 10);
```

## API Reference

| Type | Description |
|---|---|
| `CorrelationResult { coefficient, lag, significance }` | Correlation at a specific lag |
| `CrossCorrelation { series_a, series_b, correlations }` | Full cross-correlation across all lags |
| `DependencyEdge { from, to, strength, lag }` | Directed edge in dependency graph |
| `DependencyGraph { edges, nodes }` | Graph of correlated streams. `from_correlations()`, `strongest_paths()` |

### Functions

| Function | Description |
|---|---|
| `cross_correlation(a, b, max_lag)` | Compute correlation at lags 0..max_lag |
| `pearson_correlation(a, b)` | Standard Pearson r |
| `DependencyGraph::from_correlations(series, threshold)` | Build graph from named series map |

## Testing

26 tests: Pearson correlation (perfect/inverse/zero), cross-correlation with lag detection, dependency graph construction, threshold filtering, path finding, edge cases.

## License

Apache-2.0
