namespace uno_palyground.PySDRGuide.Pipeline;

public enum TvSystem { PAL_I, PAL_DK }

public readonly record struct SystemProfile(
    double LumaCutoffHz,
    double ChromaLowHz,
    double ChromaHighHz,
    double AudioOffsetHz,
    double VideoBandwidthHz
)
{
    public static SystemProfile For(TvSystem s) => s switch
    {
        TvSystem.PAL_DK => new(
            LumaCutoffHz: 5.0e6,
            ChromaLowHz: 3.9e6,
            ChromaHighHz: 4.95e6,
            AudioOffsetHz: 6.5e6,
            VideoBandwidthHz: 6.0e6
        ),
        TvSystem.PAL_I => new(
            LumaCutoffHz: 4.8e6,
            ChromaLowHz: 3.9e6,
            ChromaHighHz: 4.95e6,
            AudioOffsetHz: 6.0e6,
            VideoBandwidthHz: 5.5e6
        ),
        _ => throw new ArgumentOutOfRangeException(nameof(s))
    };
}

public enum FieldOrder
{
    TopFieldFirst,
    BottomFieldFirst
}
