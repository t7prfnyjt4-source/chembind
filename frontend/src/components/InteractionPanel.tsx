interface Interaction {
  residue: string;
  type: string;
}

const TYPE_COLORS: Record<string, string> = {
  "HBDonor": "#3b82f6",
  "HBAcceptor": "#3b82f6",
  "Hydrophobic": "#f59e0b",
  "PiStacking": "#8b5cf6",
  "PiCation": "#ec4899",
  "SaltBridge": "#ef4444",
  "VdWContact": "#6b7280",
};

export default function InteractionPanel({ interactions }: { interactions: Interaction[] }) {
  if (!interactions || interactions.length === 0) {
    return <div style={{ fontSize: 12, opacity: 0.6 }}>None detected</div>;
  }

  // Group by type
  const byType: Record<string, string[]> = {};
  for (const int of interactions) {
    if (!byType[int.type]) byType[int.type] = [];
    byType[int.type].push(int.residue);
  }

  return (
    <div style={{ fontSize: 13 }}>
      <div style={{ fontWeight: 600, marginBottom: 4 }}>Interactions</div>

      <div style={{ display: "flex", flexWrap: "wrap", gap: 4, marginBottom: 8 }}>
        {Object.entries(byType).map(([type, residues]) => (
          <span
            key={type}
            style={{
              padding: "2px 8px",
              borderRadius: 12,
              background: TYPE_COLORS[type] || "#6b7280",
              color: "#fff",
              fontSize: 11,
              fontWeight: 600,
            }}
          >
            {type} ×{residues.length}
          </span>
        ))}
      </div>

      <div style={{ fontSize: 12 }}>
        {interactions.map((int, i) => (
          <div key={i} style={{ display: "flex", gap: 8, padding: "1px 0" }}>
            <span style={{ fontFamily: "monospace", minWidth: 80 }}>{int.residue}</span>
            <span style={{ opacity: 0.7 }}>{int.type}</span>
          </div>
        ))}
      </div>
    </div>
  );
}
