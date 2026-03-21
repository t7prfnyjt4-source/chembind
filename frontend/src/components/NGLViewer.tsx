import { useEffect, useRef } from "react";
import * as NGL from "ngl";

/**
 * Convert atoms array to SDF string for NGL loading.
 */
export function sdfFromAtoms(
  atoms: { symbol: string; x: number; y: number; z: number }[]
): string {
  const n = atoms.length;
  const lines = [
    "", "ChemBind", "",
    `${String(n).padStart(3)}${String(0).padStart(3)}  0  0  0  0  0  0  0  0999 V2000`,
  ];

  for (const a of atoms) {
    lines.push(
      `${a.x.toFixed(4).padStart(10)}${a.y.toFixed(4).padStart(10)}${a.z.toFixed(4).padStart(10)} ${a.symbol.padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`
    );
  }

  lines.push("M  END", "$$$$", "");
  return lines.join("\n");
}

export default function NGLViewer({
  sdfData,
  width = 400,
  height = 400,
}: {
  sdfData: string;
  width?: number;
  height?: number;
}) {
  const containerRef = useRef<HTMLDivElement>(null);
  const stageRef = useRef<any>(null);

  useEffect(() => {
    if (!containerRef.current) return;

    const stage = new NGL.Stage(containerRef.current, {
      backgroundColor: "white",
    });
    stageRef.current = stage;

    const blob = new Blob([sdfData], { type: "text/plain" });
    stage.loadFile(blob, { ext: "sdf" }).then((comp: any) => {
      comp.addRepresentation("ball+stick");
      comp.autoView();
    });

    return () => {
      stage.dispose();
      stageRef.current = null;
    };
  }, [sdfData]);

  return (
    <div
      ref={containerRef}
      style={{
        width,
        height,
        border: "1px solid #555",
        borderRadius: 4,
      }}
    />
  );
}
