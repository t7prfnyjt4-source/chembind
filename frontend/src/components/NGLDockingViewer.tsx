import { useEffect, useRef, useState } from "react";
import * as NGL from "ngl";

interface PoseAtom {
  name: string;
  x: number;
  y: number;
  z: number;
  element: string;
}

interface Pose {
  pose_id: number;
  score: number;
  atoms: PoseAtom[];
  interactions?: { residue: string; type: string }[];
}

function atomsToSDF(atoms: PoseAtom[]): string {
  const n = atoms.length;
  const lines = [
    "", "Docking Pose", "",
    `${String(n).padStart(3)}${String(0).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`,
  ];
  for (const a of atoms) {
    const el = a.element || a.name[0];
    lines.push(
      `${a.x.toFixed(4).padStart(10)}${a.y.toFixed(4).padStart(10)}${a.z.toFixed(4).padStart(10)} ${el.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`
    );
  }
  lines.push("M  END", "$$$$", "");
  return lines.join("\n");
}

export default function NGLDockingViewer({
  poses,
  proteinPdb,
}: {
  poses: Pose[];
  proteinPdb?: string;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);
  const [selectedPose, setSelectedPose] = useState(0);

  useEffect(() => {
    if (!containerRef.current) return;

    const stage = new NGL.Stage(containerRef.current, {
      backgroundColor: "white",
    });
    stageRef.current = stage;

    // Load protein if available
    if (proteinPdb) {
      const pdbBlob = new Blob([proteinPdb], { type: "text/plain" });
      stage.loadFile(pdbBlob, { ext: "pdb" }).then((comp: any) => {
        comp.addRepresentation("cartoon", { color: "gray", opacity: 0.7 });
        comp.autoView();
      });
    }

    return () => {
      stage.dispose();
      stageRef.current = null;
    };
  }, [proteinPdb]);

  // Update ligand on pose change
  useEffect(() => {
    if (!stageRef.current || !poses.length) return;
    const stage = stageRef.current;

    // Remove ligand components (keep protein as first)
    const components = stage.compList;
    while (components.length > 1) {
      stage.removeComponent(components[components.length - 1]);
    }

    const sdf = atomsToSDF(poses[selectedPose].atoms);
    const blob = new Blob([sdf], { type: "text/plain" });
    stage.loadFile(blob, { ext: "sdf" }).then((comp: any) => {
      comp.addRepresentation("ball+stick", { colorScheme: "element" });
      if (!proteinPdb) comp.autoView();
    });
  }, [selectedPose, poses, proteinPdb]);

  const current = poses[selectedPose];

  return (
    <div>
      <div
        ref={containerRef}
        style={{ width: 500, height: 400, border: "1px solid #555", borderRadius: 4 }}
      />

      <div style={{ marginTop: 8, display: "flex", gap: 8, alignItems: "center", fontSize: 13 }}>
        <span>Pose:</span>
        {poses.map((p, i) => (
          <button
            key={p.pose_id}
            onClick={() => setSelectedPose(i)}
            style={{
              fontWeight: i === selectedPose ? 700 : 400,
              fontSize: 12,
            }}
          >
            #{p.pose_id} ({p.score.toFixed(1)})
          </button>
        ))}
      </div>

      {current?.interactions && current.interactions.length > 0 && (
        <div style={{ marginTop: 4, fontSize: 12 }}>
          Interactions: {current.interactions.map((int, i) => (
            <span key={i} style={{ padding: "1px 4px", margin: 1, background: "#e0e7ff", borderRadius: 3 }}>
              {int.residue}:{int.type}
            </span>
          ))}
        </div>
      )}
    </div>
  );
}
