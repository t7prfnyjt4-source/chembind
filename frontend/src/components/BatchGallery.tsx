import { useState, useEffect, useRef, useCallback } from "react";

interface BatchItem {
  smiles: string;
  ok: boolean;
  descriptors?: { mw?: number; logp?: number };
  thumbnail?: string;
}

type SortKey = "mw" | "logp" | "smiles";
type FilterMode = "all" | "valid";

export default function BatchGallery({
  items,
  onSelectSmiles,
}: {
  items: BatchItem[];
  onSelectSmiles?: (smiles: string) => void;
}) {
  const [filter, setFilter] = useState<FilterMode>("all");
  const [sortBy, setSortBy] = useState<SortKey>("mw");
  const [visibleCount, setVisibleCount] = useState(40);
  const sentinelRef = useRef<HTMLDivElement>(null);

  const filtered = items.filter((item) => {
    if (filter === "valid") return item.ok;
    return true;
  });

  const sorted = [...filtered].sort((a, b) => {
    if (sortBy === "smiles") return (a.smiles || "").localeCompare(b.smiles || "");
    const va = a.descriptors?.[sortBy] ?? 0;
    const vb = b.descriptors?.[sortBy] ?? 0;
    return (va as number) - (vb as number);
  });

  const visible = sorted.slice(0, visibleCount);

  // IntersectionObserver for virtual scroll
  const handleObserver = useCallback(
    (entries: IntersectionObserverEntry[]) => {
      if (entries[0].isIntersecting && visibleCount < sorted.length) {
        setVisibleCount((c) => Math.min(c + 40, sorted.length));
      }
    },
    [visibleCount, sorted.length]
  );

  useEffect(() => {
    const el = sentinelRef.current;
    if (!el) return;
    const observer = new IntersectionObserver(handleObserver, { threshold: 0.1 });
    observer.observe(el);
    return () => observer.disconnect();
  }, [handleObserver]);

  return (
    <div>
      <div style={{ display: "flex", gap: 8, marginBottom: 12, alignItems: "center", fontSize: 13 }}>
        <label>Filter:</label>
        <select value={filter} onChange={(e) => setFilter(e.target.value as FilterMode)}>
          <option value="all">All</option>
          <option value="valid">Valid only</option>
        </select>
        <label style={{ marginLeft: 12 }}>Sort:</label>
        <select value={sortBy} onChange={(e) => setSortBy(e.target.value as SortKey)}>
          <option value="mw">MW</option>
          <option value="logp">LogP</option>
          <option value="smiles">SMILES</option>
        </select>
        <span style={{ marginLeft: "auto", opacity: 0.7 }}>
          {filtered.length} items
        </span>
      </div>

      <div
        style={{
          display: "grid",
          gridTemplateColumns: "repeat(auto-fill, minmax(180px, 1fr))",
          gap: 8,
        }}
      >
        {visible.map((item, i) => (
          <div
            key={i}
            onClick={() => onSelectSmiles?.(item.smiles)}
            style={{
              border: "1px solid #ddd",
              borderRadius: 6,
              padding: 8,
              cursor: onSelectSmiles ? "pointer" : "default",
              background: item.ok ? "#fff" : "#fff5f5",
            }}
          >
            {item.thumbnail && (
              <img
                src={`data:image/png;base64,${item.thumbnail}`}
                alt={item.smiles}
                style={{ width: "100%", height: 120, objectFit: "contain" }}
              />
            )}
            <div style={{ fontFamily: "monospace", fontSize: 11, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}>
              {item.smiles}
            </div>
            <div style={{ display: "flex", gap: 6, fontSize: 11, opacity: 0.7, marginTop: 2 }}>
              <span style={{
                padding: "1px 4px", borderRadius: 3,
                background: item.ok ? "#dcfce7" : "#fee2e2",
                color: item.ok ? "#166534" : "#991b1b",
              }}>
                {item.ok ? "valid" : "invalid"}
              </span>
              {item.descriptors?.mw !== undefined && (
                <span>MW: {item.descriptors.mw.toFixed(1)}</span>
              )}
            </div>
          </div>
        ))}
      </div>

      <div ref={sentinelRef} style={{ height: 1 }} />
    </div>
  );
}
